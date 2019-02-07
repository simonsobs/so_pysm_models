import numpy as np
import healpy as hp

try:
    from pixell import curvedsky, enmap
except:
    pass
import pysm


class PrecomputedAlms:
    def __init__(
        self,
        filename,
        input_units="uK_RJ",
        input_reference_frequency_GHz=None,
        target_nside=None,
        target_shape=None,
        target_wcs=None,
        precompute_output_map=True,
        has_polarization=True,
        pixel_indices=None,
    ):
        """Generic component based on Precomputed Alms

        Load a set of Alms from a FITS file and generate maps at the requested
        resolution and frequency assuming the CMB black body spectrum.
        A single set of Alms is used for all frequencies requested by PySM,
        consider that PySM expects the output of components to be in uK_RJ.

        See more details at https://so-pysm-models.readthedocs.io/en/latest/so_pysm_models/models.html

        Parameters
        ----------
        filename : string
            Path to the input Alms in FITS format
        input_units : string
            Input unit strings as defined by pysm.convert_units, e.g. K_CMB, uK_RJ, MJysr
        input_reference_frequency_GHz : float
            If input units are K_RJ or Jysr, the reference frequency
        target_nside : int
            HEALPix NSIDE of the output maps
        precompute_output_map : bool
            If True (default), Alms are transformed into a map in the constructor,
            if False, the object only stores the Alms and generate the map at each
            call of the signal method, this is useful to generate maps convolved
            with different beams
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
        if not input_units.endswith("CMB") and input_reference_frequency_GHz is None:
            raise Exception(
                "If the input maps are in not in K_CMB, you need to specify `input_reference_frequency_GHz`"
            )
        self.input_reference_frequency_GHz = input_reference_frequency_GHz
        self.pixel_indices = pixel_indices
        self.has_polarization = has_polarization

        alm = np.complex128(
            hp.read_alm(self.filename, hdu=(1, 2, 3) if self.has_polarization else 1)
        )

        if precompute_output_map:
            self.output_map = self.compute_output_map(alm)
        else:
            self.alm = alm

    def compute_output_map(self, alm):

        if self.nside is None:
            assert (self.shape is not None) and (self.wcs is not None)
            n_comp = 3 if self.has_polarization else 1
            output_map = enmap.empty((n_comp,) + self.shape[-2:], self.wcs)
            curvedsky.alm2map(alm, self.output_map, spin=[0, 2], verbose=True)
        elif self.nside is not None:
            output_map = hp.alm2map(alm, self.nside)
        else:
            raise ValueError("You must specify either nside or both of shape and wcs")
        return output_map

    def signal(self, nu=[148.], fwhm_arcmin=None, output_units="uK_RJ", **kwargs):
        """Return map in uK_RJ at given frequency or array of frequencies

        If nothing is specified for nu, we default to providing an unmodulated map
        at 148 GHz. The value 148 Ghz does not matter if the output is in
        uK_RJ.

        Parameters
        ----------
        nu : list or ndarray
            Frequency or frequencies in GHz at which compute the signal
        fwhm_arcmin : float (optional)
            Smooth the input alms before computing the signal, this can only be used
            if the class was initialized with `precompute_output_map` to False.
        output_units : str
            Output units, as defined in `pysm.convert_units`, by default this is
            "uK_RJ" as expected by PySM.

        Returns
        -------
        output_maps : ndarray
            Output maps array with the shape (num_freqs, 1 or 3 (I or IQU), npix)
        """

        try:
            nnu = len(nu)
        except TypeError:
            nnu = 1
            nu = np.array([nu])

        try:
            output_map = self.output_map
        except AttributeError:
            if fwhm_arcmin is None:
                alm = self.alm
            else:
                alm = hp.smoothalm(
                    self.alm, fwhm=np.radians(fwhm_arcmin / 60), pol=True, inplace=False
                )

            output_map = self.compute_output_map(alm)

        # use tile to output the same map for all frequencies
        out = np.tile(output_map, (nnu, 1, 1))
        if self.wcs is not None:
            out = enmap.enmap(out, self.wcs)
        out *= (
            (
                pysm.convert_units(
                    self.input_units, "uK_CMB", self.input_reference_frequency_GHz
                )
                * pysm.convert_units("uK_CMB", output_units, nu)
            )
            .reshape((nnu, 1, 1))
            .astype(float)
        )

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out
