import os
import numpy as np
from scipy.interpolate import interp1d

import healpy as hp

import pysm


class InterpolatingComponent:

    def __init__(
        self,
        path,
        input_units,
        target_nside,
        interpolation_kind="linear",
        has_polarization=True,
        pixel_indices=None,
        mpi_comm=None,
        verbose=False,
    ):
        """PySM component interpolating between precomputed maps

        See more details at https://so-pysm-models.readthedocs.io/en/latest/so_pysm_models/models.html

        Parameters
        ----------
        path : str
            Path should contain maps named as the frequency in GHz e.g. 20.fits or 20.5.fits or 00100.fits
        input_units : str
            Any unit available in PySM (see `pysm.convert_units` e.g. `Jysr`, `MJsr`, `uK_RJ`, `K_CMB`).
        target_nside : int
            HEALPix NSIDE of the output maps
        has_polarization : bool
            whether or not to simulate also polarization maps
        pixel_indices : ndarray of ints
            Outputs partial maps given HEALPix pixel indices in RING ordering
        mpi_comm : mpi4py communicator
            See the documentation of pysm.read_map
        verbose : bool
            Control amount of output
        """

        self.maps = {}
        self.maps = self.get_filenames(path)

        self.freqs = np.array(list(self.maps.keys()))
        self.freqs.sort()
        self.input_units = input_units
        self.nside = target_nside
        self.pixel_indices = pixel_indices
        self.has_polarization = has_polarization
        self.interpolation_kind = interpolation_kind
        self.mpi_comm = mpi_comm
        self.verbose = verbose

    def get_filenames(self, path):
        # Override this to implement name convention
        filenames = {}
        for f in os.listdir(path):
            if f.endswith(".fits"):
                freq = float(os.path.splitext(f)[0])
                filenames[freq] = os.path.join(path, f)
        return filenames
        
    def signal(self, nu, **kwargs):
        """Return map at given frequency or array of frequencies"""

        if np.isscalar(nu):

            # special case: we request only 1 frequency and that is among the ones
            # available as input
            check_isclose = np.isclose(self.freqs, nu)
            if np.any(check_isclose):

                freq = self.freqs[check_isclose][0]
                out = self.read_map(freq)
                if self.has_polarization:
                    return out
                else:
                    zeros = np.zeros_like(out)
                    return np.array([out, zeros, zeros])

            else:  # continue with interpolation as with an array of nus
                nu = np.array([nu])
        else:
            nu = np.asarray(nu)

        assert (
            nu[0] >= self.freqs[0]
        ), "Frequency not supported, requested {} Ghz < lower bound {} GHz".format(
            nu[0], self.freqs[0]
        )
        assert (
            nu[-1] <= self.freqs[-1]
        ), "Frequency not supported, requested {} Ghz > upper bound {} GHz".format(
            nu[-1], self.freqs[-1]
        )

        first_freq_i, last_freq_i = np.searchsorted(self.freqs, [nu[0], nu[-1]])
        first_freq_i -= 1
        last_freq_i += 1

        freq_range = self.freqs[first_freq_i:last_freq_i]

        if self.interpolation_kind == 'linear':
            # only select frequencies that are actually necessary if linear
            nuin=self.freqs
            usenu = np.zeros(len(nuin),dtype=bool)
            for iout in np.arange(len(nu)):
                for iin in np.arange(len(nuin)-1):
                    if nuin[iin] < nu[iout] <= nuin[iin+1]:
                        usenu[iin]=usenu[iin+1]=True
                        
            freq_range = nuin[usenu]

        if self.verbose:
            print("Frequencies considered:", freq_range)

        npix = (
            len(self.pixel_indices)
            if self.pixel_indices is not None
            else hp.nside2npix(self.nside)
        )

        # allocate a single array for all maps to be used by the interpolator
        # always use size 3 for polarization because PySM always expects IQU maps

        all_maps = np.zeros(
            (len(freq_range), 3 if self.has_polarization else 1, npix), dtype=np.double
        )

        for i, freq in enumerate(freq_range):
            if self.has_polarization:
                all_maps[i] = self.read_map(freq)
                if self.verbose:
                    for i_pol, pol in enumerate("IQU"):
                        print(
                            "Mean emission at {} GHz in {}: {:.4g} uK_RJ".format(
                                freq, pol, all_maps[i][i_pol].mean()
                            )
                        )
            else:
                all_maps[i][0] = self.read_map(freq)

        out = interp1d(freq_range, all_maps, axis=0, kind=self.interpolation_kind)(nu)

        # the output of out is always 3D, (num_freqs, IQU, npix), if num_freqs is one
        # we return only a 2D array.
        if len(out) == 1:
            return out[0]
        else:
            return out

    def read_map(self, freq):
        filename = self.maps[freq]
        return self.read_map_file(freq, filename)

    def read_map_file(self, freq, filename):
        if self.verbose:
            print("Reading map {}".format(filename))
        m = pysm.read_map(
            filename,
            nside=self.nside,
            field=(0, 1, 2) if self.has_polarization else 0,
            pixel_indices=self.pixel_indices,
            mpi_comm=self.mpi_comm,
        )
        return m * pysm.convert_units(self.input_units, "uK_RJ", freq)
