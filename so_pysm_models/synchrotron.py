import numpy as np

import healpy as hp

from . import laws


class GaussianSynchrotron:

    def __init__(
        self,
        target_nside,
        has_polarization=True,
        pixel_indices=None,
        EE_amplitude=4.3,
        EtoB=0.5,
        alpha=-1.0,
        beta=-3.1,
        curv=0.,
        nu_0=23.,
        verbose=False,
    ):
        """Gaussian synchrotron model

        Parameters
        ----------
        target_nside : int
            HEALPix NSIDE of the output maps
        pixel_indices : ndarray of ints
            Outputa partial maps given HEALPix pixel indices in RING ordering
        EE_amplitude : float
            amplitude of the synchrotron EE power spectrum (D_ell) at the reference
            frequency and ell=80, in muK^2 and thermodinamic units.
            Default: 4.3 from the amplitude of S-PASS E-modes power spectrum at 2.3GHz
            in the region covered by SO-SAT, rescaled at 23GHz with a powerlaw with
            beta_s = -3.1
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
        """

        self.target_nside = target_nside
        self.pixel_indices = pixel_indices
        self.has_polarization = has_polarization
        self.EE_amplitude = EE_amplitude
        self.EtoB = EtoB
        self.alpha = alpha
        self.beta = beta
        self.curv = curv
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
        clTT_sync = clZERO
        clTE_sync = clZERO
        clEE_sync = (
            dl_prefac
            * self.EE_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        BB_amplitude = self.EE_amplitude * self.EtoB
        clBB_sync = (
            dl_prefac
            * BB_amplitude
            * ((ell + 0.1) / 80.) ** self.alpha
            * laws.black_body_cmb(self.nu_0) ** 2
        )
        spec_sync = laws.curved_power_law(nu, self.nu_0, self.beta, self.curv)
        clt_sync = np.zeros([3, nnu, nnu, nell])
        clt_sync[0] = (
            clTT_sync[None, None, :]
            * spec_sync[:, None, None]
            * spec_sync[None, :, None]
        )
        clt_sync[1] = (
            clEE_sync[None, None, :]
            * spec_sync[:, None, None]
            * spec_sync[None, :, None]
        )
        clt_sync[2] = (
            clBB_sync[None, None, :]
            * spec_sync[:, None, None]
            * spec_sync[None, :, None]
        )
        # Generate foreground maps
        amp_sync = np.array(
            hp.synfast(
                [clTT_sync, clEE_sync, clBB_sync, clTE_sync, clZERO, clZERO],
                self.target_nside,
                pol=True,
                new=True,
                verbose=False,
            )
        )
        out = amp_sync[None, :, :] * spec_sync[:, None, None]

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out
