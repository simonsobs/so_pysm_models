import numpy as np

import healpy as hp

import math

from . import laws
from . import filter_utils

try:  # PySM >= 3.2.1
    import pysm3.units as u
    import pysm3 as pysm
except ImportError:
    import pysm.units as u
    import pysm


class GaussianSynchrotron(pysm.Model):
    def __init__(
        self,
        nside,
        has_polarization=True,
        pixel_indices=None,
        TT_amplitude=20.0,
        Toffset=72.0,
        EE_amplitude=4.3,
        rTE=0.35,
        EtoB=0.5,
        alpha=-1.0,
        beta=-3.1,
        curv=0.0,
        nu_0=23.0,
        seed=None,
        map_dist=None,
    ):
        """Gaussian synchrotron model

        See more details at https://so-pysm-models.readthedocs.io/en/latest/so_pysm_models/models.html

        Parameters
        ----------
        nside : int
            HEALPix NSIDE of the output maps
        has_polarization : bool
            whether or not to simulate also polarization maps
            Default: True
        pixel_indices : ndarray of ints
            Outputa partial maps given HEALPix pixel indices in RING ordering
        TT_amplitude : float
            amplitude of synchrotron TT power spectrum (D_ell) at at the reference
            frequency and ell=80, in muK^2 and thermodinamic units.
            Default: 20 from the amplitude of PySM-s0 synchrotron model at 23GHz
            in the region covered by SO-SAT.
        Toffset : float
            offset to be applied to the temperature map in muK.
            Default: 72 from the mean value of the T PySM-s0 synch map at 23GHz
            in the region covered by SO-SAT
        EE_amplitude : float
            same as TT_amplitude but for EE power spectrum.
            Default: 4.3 from the amplitude of S-PASS E-modes power spectrum at 2.3GHz
            in the region covered by SO-SAT, rescaled at 23GHz with a powerlaw with
            beta_s = -3.1
        rTE : float
            TE correlation factor defined as: rTE = clTE/sqrt(clTT*clEE)
            Default: 0.35 from Planck IX 2018
        EtoB : float
            ratio between E and B-mode amplitude.
            Default: 0.5 from Krachmalnicoff et al. 2018
        alpha : spectral tilt of the synchrotron power spectrum (D_ell).
            Default: -1.0 from Krachmalnicoff et al. 2018
        beta : synchrotron spectral index.
            Default: -3.1 from Planck 2018 IX
        curv : synchrotron curvature index.
            Default: 0.
        nu_0 : synchrotron reference frequency in GHz.
            Default: 23
        seed : int
            seed for random realization of map
            Default: None
        """

        super().__init__(nside=nside, map_dist=map_dist)
        self.has_polarization = has_polarization
        self.TT_amplitude = TT_amplitude
        self.Toffset = Toffset
        self.EE_amplitude = EE_amplitude
        self.EtoB = EtoB
        self.alpha = alpha
        self.beta = beta
        self.curv = curv
        self.nu_0 = nu_0
        self.seed = seed

    @u.quantity_input
    def get_emission(self, freqs: u.GHz, weights=None) -> u.uK_RJ:
        """Return map in uK_RJ at given frequency or array of frequencies"""

        nu = pysm.check_freq_input(freqs)
        assert (
            len(nu) == 1
        ), "Bandpass integration not implemented in Gaussian emissions"
        nell = 3 * self.nside
        ell = np.arange(nell)
        dl_prefac = 2 * np.pi / ((ell + 0.01) * (ell + 1))
        clZERO = np.zeros_like(ell)
        clTT_sync = (
            dl_prefac
            * self.TT_amplitude
            * ((ell + 0.1) / 80.0) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        clTT_sync[0] = 0.0
        clEE_sync = (
            dl_prefac
            * self.EE_amplitude
            * ((ell + 0.1) / 80.0) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        BB_amplitude = self.EE_amplitude * self.EtoB
        clBB_sync = (
            dl_prefac
            * BB_amplitude
            * ((ell + 0.1) / 80.0) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        if self.seed == None:
            mseed = np.random.randint(0, 2 ** 32)
        else:
            mseed = self.seed
        np.random.seed(mseed)
        if self.nside <= 64:
            nside_temp = self.nside
        else:
            nside_temp = 64
        amp_sync = np.array(
            hp.synfast(
                [clTT_sync, clEE_sync, clBB_sync, clZERO, clZERO, clZERO],
                nside_temp,
                pol=True,
                new=True,
                verbose=False,
            )
        )
        amp_sync[0] += self.Toffset
        lbreak_TT = 2
        while np.any(amp_sync[0] < 0):
            clTT_sync[1:lbreak_TT] = clTT_sync[lbreak_TT]
            np.random.seed(mseed)
            amp_sync = np.array(
                hp.synfast(
                    [clTT_sync, clEE_sync, clBB_sync, clZERO, clZERO, clZERO],
                    nside_temp,
                    pol=True,
                    new=True,
                    verbose=False,
                )
            )
            amp_sync[0] += self.Toffset
            lbreak_TT += 1
        if self.nside > 64:
            low_pass_filter = filter_utils.create_low_pass_filter(
                l1=30, l2=60, lmax=64 * 3 - 1
            )
            high_pass_filter = filter_utils.create_high_pass_filter(
                l1=30, l2=60, lmax=nell - 1
            )
            clTT_sync_hpf = clTT_sync * high_pass_filter
            amp_sync[0] = filter_utils.apply_filter(amp_sync[0], low_pass_filter)
            np.random.seed(mseed)
            amp_sync_hell = np.array(
                hp.synfast(
                    [clTT_sync_hpf, clEE_sync, clBB_sync, clZERO, clZERO, clZERO],
                    self.nside,
                    pol=True,
                    new=True,
                    verbose=False,
                )
            )
            amp_sync = hp.ud_grade(amp_sync, self.nside)
            amp_sync[1:3] = amp_sync[1:3] * 0.0
            amp_sync += amp_sync_hell
        min_map = np.min(amp_sync[0])
        if min_map < 0:
            Toffset_add = math.ceil(-min_map)
            amp_sync[0] += Toffset_add
        spec_sync = laws.curved_power_law(nu, self.nu_0, self.beta, self.curv)
        out = amp_sync[None, :, :] * spec_sync[:, None, None]
        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        out = out[0]
        return out << u.uK_RJ
