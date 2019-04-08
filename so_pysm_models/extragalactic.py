from . import InterpolatingComponent
from . import utils
import pysm
import numpy as np
import healpy as hp


def y2uK_CMB(nu):
    """Compton-y distortion at a given frequency

    Parmeters:
    nu (float): frequency in GHz

    Returns:
    float: intensity variation dT_CMB in micro-Kelvin
      dT_CMB = dI_nu / (dB_nu / dT)_Tcmb
      where B_nu is the Planck function and dI_nu is the intensity distortion

    """

    h = 6.62607004e-27
    k = 1.380622e-16
    Tcmb = 2.725
    x = h * nu * 1e9 / k / Tcmb
    return 1e6 * Tcmb * (x * (np.exp(x) + 1) / (np.exp(x) - 1) - 4)


class WebSkyCIB(InterpolatingComponent):
    """PySM component interpolating between precomputed maps"""

    def __init__(
        self,
        websky_version="0.3",
        input_units="MJysr",
        target_nside=4096,
        interpolation_kind="linear",
        pixel_indices=None,
        mpi_comm=None,
        verbose=False,
    ):
        super().__init__(
            websky_version,
            input_units,
            target_nside,
            interpolation_kind,
            has_polarization=False,
            pixel_indices=pixel_indices,
            mpi_comm=mpi_comm,
            verbose=verbose,
        )

    def get_filenames(self, path):
        """Get filenames for a websky version
        For a standard interpolating component, we list files in folder,
        here we need to know the names in advance so that we can only download the required maps
        """

        websky_version = path
        if websky_version == "0.3":

            available_frequencies = []
            for base_freq in [27, 39, 93, 145, 225, 280]:
                for delta_freq in [-1, 0, 1]:
                    available_frequencies.append(base_freq + delta_freq)
            for base_freq in [100, 217, 353, 545, 857]:
                delta_freq = 0
                available_frequencies.append(base_freq + delta_freq)

            filenames = {
                freq: "websky/0.3/cib_{:04d}.fits".format(freq)
                for freq in available_frequencies
            }

        return filenames

    def read_map(self, freq):
        filename = utils.get_data_from_url(self.maps[freq])
        return self.read_map_file(freq, filename)


class WebSkySZ:

    def __init__(
        self,
        version="0.3",
        sz_type="kinetic",
        target_nside=4096,
        pixel_indices=None,
        mpi_comm=None,
        verbose=False,
    ):

        self.version = version
        self.sz_type = sz_type
        self.nside = target_nside
        self.pixel_indices = pixel_indices
        self.mpi_comm = mpi_comm
        self.verbose = verbose

    def get_filename(self):
        """Get SZ filenames for a websky version"""

        version = self.version

        if self.sz_type == "kinetic":
            filename = "websky/" + version + "/ksz.fits"
        elif self.sz_type == "thermal":
            filename = "websky/" + version + "/tsz.fits"

        return filename

    def signal(self, nu, **kwargs):
        """Return map in uK_RJ at given frequency or array of frequencies"""

        if np.isscalar(nu):
            nu = np.array([nu])

        filename = utils.get_data_from_url(self.get_filename())
        m = pysm.read_map(
            filename,
            nside=self.nside,
            field=0,
            pixel_indices=self.pixel_indices,
            mpi_comm=self.mpi_comm,
        )

        npix = (
            len(self.pixel_indices)
            if self.pixel_indices is not None
            else hp.nside2npix(self.nside)
        )

        all_maps = np.zeros((len(nu), 1, npix), dtype=np.double)

        szfac = np.ones(len(nu))
        if self.sz_type == "thermal":
            szfac = y2uK_CMB(nu)

        all_maps[:, 0, :] = np.outer(
            pysm.convert_units("uK_CMB", "uK_RJ", nu) * szfac, m
        )

        return all_maps
