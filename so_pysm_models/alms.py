import numpy as np
import healpy as hp
try: from pixell import curvedsky
except: pass

import pysm

class PrecomputedAlms:

    def __init__(
        self, filename, input_units,
            target_nside = None, 
            target_shape = None,
            target_wcs = None,
            has_polarization=True, pixel_indices=None
    ):
        """Generic component based on Precomputed Alms

        A single set of Alms is used for all frequencies requested by PySM,
        consider that PySM expects the output of components to be in uK_RJ.

        Parameters
        ----------
        target_nside : int
            HEALPix NSIDE of the output maps
        filename : string
            Path to the input Alms in FITS format
        input_units : string
            Input unit strings as defined by pysm.convert_units, e.g. K_CMB, uK_RJ, MJysr
        has_polarization : bool
            whether or not to simulate also polarization maps
            Default: True
        pixel_indices : ndarray of ints
            Output a partial maps given HEALPix pixel indices in RING ordering
        """

        self.nside = target_nside
        self.shape = target_shape
        self.wcs = target_wcs
        self.filename = filename
        self.input_units = input_units
        self.pixel_indices = pixel_indices
        self.has_polarization = has_polarization

        self.alm = np.complex128(
            hp.read_alm(self.filename, hdu=(1, 2, 3) if self.has_polarization else 1)
        )

        # use tile to output the same map for all frequencies
        if self.nside is None:
            ncomp = 3 if self.has_polarization else 1
            self.omap = enmap.empty( (ncomp,)+self.shape[-2:], self.wcs)
            curvedsky.alm2map(self.alm, self.omap, spin = [0, 2], verbose = True)
        else:
            self.omap = hp.alm2map(self.alm, self.nside)
        

    def signal(self, nu, **kwargs):
        """Return map in uK_RJ at given frequency or array of frequencies"""

        try:
            nnu = len(nu)
        except TypeError:
            nnu = 1
            nu = np.array([nu])
        
        out = np.tile(self.omap, (nnu, 1, 1))
        out *= pysm.convert_units(self.input_units, "uK_RJ", nu).reshape((nnu, 1, 1))

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out
