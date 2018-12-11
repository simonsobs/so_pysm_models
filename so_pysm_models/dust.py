import numpy as np

import healpy as hp

from . import laws


class GaussianDust:

    def __init__(
        self,
        target_nside,
        has_polarization=True,
        pixel_indices=None,
        TT_amplitude=350.,
        Toffset=18.,
        EE_amplitude=100.,
        rTE=0.35,
        EtoB=0.5,
        alpha=-0.42,
        beta=1.53,
        temp=19.6,
        nu_0=353,
        seed=None,
    ):
        """Gaussian dust model

        Parameters
        ----------
        target_nside : int
            HEALPix NSIDE of the output maps
        has_polarization : bool
            whether or not to simulate also polarization maps
            Default: True
        pixel_indices : ndarray of ints
            Outputs partial maps given HEALPix pixel indices in RING ordering
        TT_amplitude : float
            amplitude of synchrotron TT power spectrum (D_ell) at at the reference
            frequency and ell=80, in muK^2 and thermodinamic units.
            Default: 350. from the amplitude of PySM-d0 dust model at 353GHz
            in the region covered by SO-SAT.
        Toffset : float
            offset to be applied to the temperature map in muK in RJ units.
            Default: 18 from the mean value of the T PySM-s0 synch map at 23GHz
            in the region covered by SO-SAT
        EE_amplitude : float
            Amplitude of EE modes D_ell at reference frequency at ell=80
            Default: 100. from the amplitude of HFI-353 E-modes spectrum in the
            region covered by SO-SAT
        rTE : float
            TE correlation factor defined as: rTE = clTE/sqrt(clTT*clEE)
            Default: 0.35 from Planck IX 2018
        EtoB: float
            ratio between E and B-mode amplitude for dust.
            Default: 0.5 from Planck 2018 IX
        alpha : same as alpha_sync for dust.
            Default: -0.42 from Planck 2018 IX
        beta : float
            dust spectral index.
            Default: 1.53 from Planck 2018 IX
        temp : float
            dust temperature.
            Default: 19.6 from Planck 2018 IX
        nu0 : float
            dust reference frequency in GHz.
            Default: 353
        seed : int
            seed for random realization of map
            Default: None
        """

        self.target_nside = target_nside
        self.pixel_indices = pixel_indices
        self.has_polarization = has_polarization
        self.TT_amplitude = TT_amplitude
        self.Toffset = Toffset
        self.EE_amplitude = EE_amplitude
        self.rTE = rTE
        self.EtoB = EtoB
        self.alpha = alpha
        self.beta = beta
        self.temp = temp
        self.nu_0 = nu_0
        self.seed = seed

    def signal(self, nu, **kwargs):
        """Return map in uK_RJ at given frequency or array of frequencies"""

        nell = 3 * self.target_nside
        try:
            nnu = len(nu)
        except TypeError:
            nnu = 1
            nu = np.array([nu])
        ell = np.arange(nell)
        dl_prefac = 2 * np.pi / ((ell + 0.01) * (ell + 1))
        clZERO = np.zeros_like(ell)
        clTT_dust = (
            dl_prefac
            * self.TT_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        clTT_dust[0] = 0.
        clEE_dust = (
            dl_prefac
            * self.EE_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        clTE_dust = self.rTE * np.sqrt(clTT_dust * clEE_dust)
        BB_amplitude = self.EE_amplitude * self.EtoB
        clBB_dust = (
            dl_prefac
            * BB_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        if self.seed == None:
            mseed = np.random.randint(0, 2 ** 32)
        else:
            mseed = self.seed
        np.random.seed(mseed)
        amp_dust = np.array(
            hp.synfast(
                [clTT_dust, clEE_dust, clBB_dust, clTE_dust, clZERO, clZERO],
                self.target_nside,
                pol=True,
                new=True,
                verbose=False,
            )
        )
        amp_dust[0] = amp_dust[0]+self.Toffset
        lbreak_TT = 2
        while np.any(amp_dust[0] < 0):
            print(lbreak_TT)
            clTT_dust[1:lbreak_TT] = clTT_dust[lbreak_TT]
            clTE_dust = self.rTE * np.sqrt(clTT_dust * clEE_dust)
            np.random.seed(mseed)
            amp_dust = np.array(
                hp.synfast(
                    [clTT_dust, clEE_dust, clBB_dust, clTE_dust, clZERO, clZERO],
                    self.target_nside,
                    pol=True,
                    new=True,
                    verbose=False,
                )
            )
            amp_dust[0] = amp_dust[0]+self.Toffset
            lbreak_TT += 1
        spec_dust = laws.modified_black_body(nu, self.nu_0, self.beta, self.temp)
        out = amp_dust[None, :, :] * spec_dust[:, None, None]

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out
