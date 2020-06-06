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


class GaussianDust(pysm.Model):
    def __init__(
        self,
        nside,
        has_polarization=True,
        pixel_indices=None,
        TT_amplitude=350.0,
        Toffset=18.0,
        EE_amplitude=100.0,
        rTE=0.35,
        EtoB=0.5,
        alpha=-0.42,
        beta=1.53,
        temp=19.6,
        nu_0=353,
        seed=None,
        map_dist=None,
    ):
        """Gaussian dust model

        See more details at https://so-pysm-models.readthedocs.io/en/latest/so_pysm_models/models.html

        Parameters
        ----------
        nside : int
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

        super().__init__(nside=nside, map_dist=map_dist)
        self.has_polarization = has_polarization
        self.TT_amplitude = TT_amplitude
        self.Toffset = Toffset
        self.EE_amplitude = EE_amplitude
        self.EtoB = EtoB
        self.alpha = alpha
        self.beta = beta
        self.temp = temp
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
        clTT_dust = (
            dl_prefac
            * self.TT_amplitude
            * ((ell + 0.1) / 80.0) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        clTT_dust[0] = 0.0
        clEE_dust = (
            dl_prefac
            * self.EE_amplitude
            * ((ell + 0.1) / 80.0) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        BB_amplitude = self.EE_amplitude * self.EtoB
        clBB_dust = (
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
        amp_dust = np.array(
            hp.synfast(
                [clTT_dust, clEE_dust, clBB_dust, clZERO, clZERO, clZERO],
                nside_temp,
                pol=True,
                new=True,
                verbose=False,
            )
        )
        amp_dust[0] = amp_dust[0] + self.Toffset
        lbreak_TT = 2
        while np.any(amp_dust[0] < 0):
            clTT_dust[1:lbreak_TT] = clTT_dust[lbreak_TT]
            np.random.seed(mseed)
            amp_dust = np.array(
                hp.synfast(
                    [clTT_dust, clEE_dust, clBB_dust, clZERO, clZERO, clZERO],
                    nside_temp,
                    pol=True,
                    new=True,
                    verbose=False,
                )
            )
            amp_dust[0] = amp_dust[0] + self.Toffset
            lbreak_TT += 1
        if self.nside > 64:
            low_pass_filter = filter_utils.create_low_pass_filter(
                l1=30, l2=60, lmax=64 * 3 - 1
            )
            high_pass_filter = filter_utils.create_high_pass_filter(
                l1=30, l2=60, lmax=nell - 1
            )
            clTT_dust_hpf = clTT_dust * high_pass_filter
            amp_dust[0] = filter_utils.apply_filter(amp_dust[0], low_pass_filter)
            np.random.seed(mseed)
            amp_dust_hell = np.array(
                hp.synfast(
                    [clTT_dust_hpf, clEE_dust, clBB_dust, clZERO, clZERO, clZERO],
                    self.nside,
                    pol=True,
                    new=True,
                    verbose=False,
                )
            )
            amp_dust = hp.ud_grade(amp_dust, self.nside)
            amp_dust[1:3] = amp_dust[1:3] * 0.0
            amp_dust += amp_dust_hell
        min_map = np.min(amp_dust[0])
        if min_map < 0:
            Toffset_add = math.ceil(-min_map)
            amp_dust[0] += Toffset_add
        spec_dust = laws.modified_black_body(nu, self.nu_0, self.beta, self.temp)
        out = amp_dust[None, :, :] * spec_dust[:, None, None]
        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        return out[0] << u.uK_RJ
