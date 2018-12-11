import numpy as np
import healpy as hp


class PrecomputedAlms:

    def __init__(
        self, target_nside, filename, has_polarization=True, pixel_indices=None
    ):
        """Generic component based on Precomputed Alms

        Parameters
        ----------
        target_nside : int
            HEALPix NSIDE of the output maps
        filename : string
            Path to the input Alms in FITS format
        has_polarization : bool
            whether or not to simulate also polarization maps
            Default: True
        pixel_indices : ndarray of ints
            Output a partial maps given HEALPix pixel indices in RING ordering
        """

        self.target_nside = target_nside
        self.filename = filename
        self.pixel_indices = pixel_indices
        self.has_polarization = has_polarization

        self.alm = np.complex128(
            hp.read_alm(self.filename, hdu=(1, 2, 3) if self.has_polarization else 1)
        )

    def signal(self, nu, **kwargs):
        """Return map in uK_RJ at given frequency or array of frequencies"""

        try:
            nnu = len(nu)
        except TypeError:
            nnu = 1
            nu = np.array([nu])

        # use tile to output the same map for all frequencies

        out = np.tile(hp.alm2map(self.alm, self.target_nside), (nnu, 1, 1))

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out
