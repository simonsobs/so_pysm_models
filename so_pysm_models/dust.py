import numpy as np

import healpy as hp

from . import laws


class GaussianDust:

    def __init__(
        self,
        target_nside,
        pixel_indices=None,
        EE_amplitude=100,
        EtoB=0.5,
        alpha=-0.42,
        beta=1.53,
        temp=19.6,
        nu_0=353,
        verbose=False,
    ):
        """Gaussian dust model

        Parameters
        ----------
        target_nside : int
            HEALPix NSIDE of the output maps
        pixel_indices : ndarray of ints
            Outputs partial maps given HEALPix pixel indices in RING ordering
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
        """

        self.target_nside = target_nside
        self.pixel_indices = pixel_indices
        self.EE_amplitude = EE_amplitude
        self.EtoB = EtoB
        self.alpha = alpha
        self.beta = beta
        self.temp = temp
        self.nu_0 = nu_0
        self.verbose = verbose

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
        clTT_dust = clZERO
        clTE_dust = clZERO
        clEE_dust = (
            dl_prefac
            * self.EE_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        BB_amplitude = self.EE_amplitude * self.EtoB
        clBB_dust = (
            dl_prefac
            * BB_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        spec_dust = laws.modified_black_body(nu, self.nu_0, self.beta, self.temp)
        clt_dust = np.zeros([3, nnu, nnu, nell])
        clt_dust[0] = (
            clTT_dust[None, None, :]
            * spec_dust[:, None, None]
            * spec_dust[None, :, None]
        )
        clt_dust[1] = (
            clEE_dust[None, None, :]
            * spec_dust[:, None, None]
            * spec_dust[None, :, None]
        )
        clt_dust[2] = (
            clBB_dust[None, None, :]
            * spec_dust[:, None, None]
            * spec_dust[None, :, None]
        )
        amp_dust = np.array(
            hp.synfast(
                [clTT_dust, clEE_dust, clBB_dust, clTE_dust, clZERO, clZERO],
                self.target_nside,
                pol=True,
                new=True,
                verbose=False,
            )
        )
        out = amp_dust[None, :, :] * spec_dust[:, None, None]

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out
