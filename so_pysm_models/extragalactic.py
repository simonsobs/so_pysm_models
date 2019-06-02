import os.path

from numba import njit
import numpy as np

from pysm import InterpolatingComponent, Model
from pysm import units as u

from pysm.utils import normalize_weights, trapz_step_inplace
from .alms import PrecomputedAlms
from . import utils



@njit
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
        input_units="MJy / sr",
        nside=4096,
        interpolation_kind="linear",
        map_dist=None,
        verbose=False,
        local_folder=None,
    ):
        self.local_folder = local_folder
        super().__init__(
            str(websky_version),
            input_units,
            nside,
            interpolation_kind,
            has_polarization=False,
            map_dist=map_dist,
            verbose=verbose,
        )
        self.dataurl = None  # utils.DATAURL

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
        if self.local_folder is not None:
            for freq in filenames:
                filenames[freq] = os.path.join(self.local_folder, filenames[freq])

        return filenames

    def read_map_by_frequency(self, freq):
        filename = utils.get_data_from_url(self.maps[freq])
        return self.read_map_file(freq, filename)


class WebSkySZ(Model):
    def __init__(
        self, version="0.3", sz_type="kinetic", nside=4096, map_dist=None, verbose=False
    ):

        super().__init__(nside=nside, map_dist=map_dist)
        self.version = str(version)
        self.sz_type = sz_type
        self.verbose = verbose

    def get_filename(self):
        """Get SZ filenames for a websky version"""

        version = self.version

        if self.sz_type == "kinetic":
            filename = "websky/" + version + "/ksz.fits"
        elif self.sz_type == "thermal":
            filename = "websky/" + version + "/tsz.fits"

        return filename

    @u.quantity_input
    def get_emission(self, freqs: u.GHz, weights=None) -> u.uK_RJ:

        nu = freqs.to(u.GHz)
        weights = normalize_weights(freqs, weights)

        if nu.isscalar:
            nu = nu.reshape(1)

        filename = utils.get_data_from_url(self.get_filename())
        m = self.read_map(filename, field=0, unit=u.uK_CMB)

        weights = (weights * u.uK_CMB).to_value(
            u.uK_RJ, equivalencies=u.cmb_equivalencies(nu)
        )

        is_thermal = self.sz_type == "thermal"
        output = (
            get_sz_emission_numba(nu.value, weights, m.value, is_thermal)
            << u.uK_RJ
        )

        # the output of out is always 2D, (IQU, npix)
        return output


@njit(parallel=True)
def get_sz_emission_numba(freqs, weights, m, is_thermal):
    output = np.zeros((3, len(m)), dtype=m.dtype)
    for i in range(len(freqs)):
        if is_thermal:
            signal = m * m.dtype.type(y2uK_CMB(freqs[i]))
        else:
            signal = m
        trapz_step_inplace(freqs, weights, i, signal, output[0])
    return output


class WebSkyCMB(PrecomputedAlms):
    def __init__(
        self,
        websky_version,
        nside,
        precompute_output_map=False,
        seed=1,
        lensed=True,
        map_dist=None,
    ):
        filename = utils.get_data_from_url(
            "websky/{}/{}lensed_alm_seed{}.fits".format(
                websky_version, "" if lensed else "un", seed
            )
        )
        super().__init__(
            filename,
            input_units="uK_CMB",
            input_reference_frequency=None,
            nside=nside,
            precompute_output_map=precompute_output_map,
            has_polarization=True,
            map_dist=map_dist,
        )
