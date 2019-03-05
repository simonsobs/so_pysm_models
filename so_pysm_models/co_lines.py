import numpy as np

import healpy as hp


class COLines:
    def __init__(
        self,
        target_nside,
        has_polarization=True,
        line="10",
        include_high_galactic_latitude_clouds=False,
        polarization_fraction=0.001,
        theta_high_galactic_latitude_deg=20.,
        random_seed=1234567,
        verbose=False,
        run_mcmole3d=False,
    ):

        """Class defining attributes for CO line emission.
           CO templates are extracted from Type 1 CO Planck maps. 
           See further details in https://www.aanda.org/articles/aa/abs/2014/11/aa21553-13/aa21553-13.html

        Parameters
        ----------
        target_nside : int
            HEALPix NSIDE of the output maps
        has_polarization : bool
            whether or not to simulate also polarization maps
        line : string 
            CO rotational transitions. 
            Accepted values : 10, 21, 32  
        polarization_fraction: float
            polarisation fraction for polarised CO emission.
        include_high_galactic_latitude_clouds: bool 
            If True it includes a simulation from MCMole3D to include 
            high Galactic Latitude clouds.
            (See more details at http://giuspugl.github.io/mcmole/index.html) 
        run_mcmole3d: bool 
            If True it simulates  HGL cluds by running MCMole3D, otherwise it coadds
            a map of HGL emission.
        random_seed: int 
            set random seed for mcmole3d simulations. 
        theta_high_galactic_latitude_deg : float
            Angle in degree  to identify High Galactic Latitude clouds 
            (i.e. clouds whose latitude b is |b|> theta_high_galactic_latitude_deg). 
        """

        self.line = line
        self.line_index = {"10": 0, "21": 1, "32": 2}[line]

        self.target_nside = target_nside

        self.template_nside = 512 if self.target_nside <= 512 else 2048
        self.planck_templatemap = hp.read_map(
            "/global/cscratch1/sd/giuspugl/CO_data/HFI_CompMap_CO-Type1_{}_R2.00_ring.fits".format(
                self.template_nside
            ),
            field=self.line_index,
            verbose=False,
        )

        self.include_high_galactic_latitude_clouds = (
            include_high_galactic_latitude_clouds
        )
        self.has_polarization = has_polarization
        self.polarization_fraction = polarization_fraction
        self.theta_high_galactic_latitude_deg = theta_high_galactic_latitude_deg
        self.random_seed = random_seed
        self.run_mcmole3d = run_mcmole3d

        self.verbose = verbose

    def signal(self):
        """
        Simulate CO signal 
        """
        out = hp.ud_grade(map_in=self.planck_templatemap, nside_out=self.target_nside)

        if self.include_high_galactic_latitude_clouds:
            out += self.simulate_high_galactic_latitude_CO()

        if self.has_polarization:
            Q_map, U_map = self.simulate_polarized_emission(out)
            return np.array([out, Q_map, U_map])
        else:
            return out

    def simulate_polarized_emission(self, I_map):
        """
        Add polarized emission by means of:
        - an overall constant polarization fraction, 
        - a depolarization map to mimick the line of sight depolarization 
          effect at low Galactic latitudes 
        - a polarization angle map coming from a dust template 
          (we exploit the observed correlation between polarized dust and 
          molecular emission in star forming regions). 
        """
        polangle = hp.read_map(
            "/global/cscratch1/sd/giuspugl/CO_data/psimap_dust90_{}.fits".format(
                self.template_nside
            ),
            verbose=False,
        )
        depolmap = hp.read_map(
            "/global/cscratch1/sd/giuspugl/CO_data/gmap_dust90_{}.fits".format(
                self.template_nside
            ),
            verbose=False,
        )

        if hp.get_nside(depolmap) != self.target_nside:
            polangle = hp.ud_grade(map_in=polangle, nside_out=self.target_nside)
            depolmap = hp.ud_grade(map_in=depolmap, nside_out=self.target_nside)

        cospolangle = np.cos(2. * polangle)
        sinpolangle = np.sin(2. * polangle)

        Q_map = self.polarization_fraction * depolmap * cospolangle * I_map
        U_map = self.polarization_fraction * depolmap * sinpolangle * I_map
        return Q_map, U_map

    def simulate_high_galactic_latitude_CO(self):
        """
        Coadd High Galactic Latitude CO emission, simulated with  MCMole3D.
        
        """
        if self.run_mcmole3d:
            import mcmole3d as cl

            # params to MCMole
            N = 40000
            L_0 = 20.4  # pc
            L_min = .3
            L_max = 60.
            R_ring = 5.8
            sigma_ring = 2.7  # kpc
            R_bulge = 3.
            R_z = 10  # kpc
            z_0 = 0.1
            Em_0 = 240.
            R_em = 6.6
            model = "LogSpiral"

            nside = self.target_nside
            Itot_o, _ = cl.integrate_intensity_map(
                self.planck_templatemap,
                hp.get_nside(self.planck_templatemap),
                planck_map=True,
            )
            Pop = cl.Cloud_Population(N, model, randseed=self.random_seed)

            Pop.set_parameters(
                radial_distr=[R_ring, sigma_ring, R_bulge],
                typical_size=L_0,
                size_range=[L_min, L_max],
                thickness_distr=[z_0, R_z],
                emissivity=[Em_0, R_em],
            )
            Pop()

            if self.verbose:
                Pop.print_parameters()
            # project into  Healpix maps
            mapclouds = cl.do_healpy_map(
                Pop,
                nside,
                highgalcut=np.deg2rad(90. - self.theta_high_galactic_latitude_deg),
                apodization="gaussian",
                verbose=self.verbose,
            )
            Itot_m, _ = cl.integrate_intensity_map(mapclouds, nside)
            # convert simulated map into the units of the Planck one
            rescaling_factor = Itot_m / Itot_o
            mapclouds /= rescaling_factor
            hglmask = np.zeros_like(mapclouds)
            # Apply mask to low galactic latitudes
            listhgl = hp.query_strip(
                nside,
                np.deg2rad(90. + self.theta_high_galactic_latitude_deg),
                np.deg2rad(90 - self.theta_high_galactic_latitude_deg),
            )
            hglmask[listhgl] = 1.
            rmsplanck = self.planck_templatemap[listhgl].std()
            rmssim = mapclouds[listhgl].std()
            if rmssim == 0.:
                belowplanck = 1.
            else:
                belowplanck = rmssim / rmsplanck

            return mapclouds * hglmask / belowplanck
        else:
            mapclouds = hp.read_map(
                "/global/cscratch1/sd/giuspugl/CO_data/mcmoleCO_HGL_{}.fits".format(
                    self.template_nside
                ),
                field=self.line_index,
                verbose=False,
            )

            if hp.get_nside(mapclouds) != self.target_nside:
                mapclouds = hp.ud_grade(map_in=mapclouds, nside_out=self.target_nside)
            return mapclouds
