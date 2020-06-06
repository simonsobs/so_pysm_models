import os.path
from pathlib import Path

from numba import njit
import numpy as np

try:  # PySM >= 3.2.1
    import pysm3.units as u
    import pysm3 as pysm
except ImportError:
    import pysm.units as u
    import pysm

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


class WebSkyCIB(pysm.InterpolatingComponent):
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
        coord="C",
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
        self.remote_data = utils.RemoteData(coord)

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
                freq: "websky/0.3/{nside}cib_{:04d}.fits".format(
                    freq, nside="512/" if self.nside <= 512 else ""
                )
                for freq in available_frequencies
            }
        if self.local_folder is not None:
            for freq in filenames:
                filenames[freq] = os.path.join(self.local_folder, filenames[freq])

        return filenames

    def read_map_by_frequency(self, freq):
        filename = self.remote_data.get(self.maps[freq])
        return self.read_map_file(freq, filename)


class WebSkySZ(pysm.Model):
    def __init__(
        self,
        version="0.3",
        sz_type="kinetic",
        nside=4096,
        map_dist=None,
        verbose=False,
        coord="C",
    ):

        super().__init__(nside=nside, map_dist=map_dist)
        self.version = str(version)
        self.sz_type = sz_type
        self.verbose = verbose
        self.remote_data = utils.RemoteData(coord)
        filename = self.remote_data.get(self.get_filename())
        self.m = self.read_map(filename, field=0, unit=u.uK_CMB)

    def get_filename(self):
        """Get SZ filenames for a websky version"""

        path = Path("websky") / self.version

        if self.nside <= 512:
            path /= "512"

        if self.sz_type == "kinetic":
            path = path / "ksz.fits"
        elif self.sz_type == "thermal":
            path = path / "tsz.fits"

        return str(path)

    @u.quantity_input
    def get_emission(self, freqs: u.GHz, weights=None) -> u.uK_RJ:

        freqs = pysm.check_freq_input(freqs)
        weights = pysm.normalize_weights(freqs, weights)

        # input map is in uK_CMB, we multiply the weights which are
        # in uK_RJ by the conversion factor of uK_CMB->uK_RJ
        # this is the equivalent of
        weights = (weights * u.uK_CMB).to_value(
            u.uK_RJ, equivalencies=u.cmb_equivalencies(freqs * u.GHz)
        )

        is_thermal = self.sz_type == "thermal"
        output = (
            get_sz_emission_numba(freqs, weights, self.m.value, is_thermal) << u.uK_RJ
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
        pysm.utils.trapz_step_inplace(freqs, weights, i, signal, output[0])
    return output


class WebSkyCMBTensor(PrecomputedAlms):
    def __init__(
        self,
        websky_version,
        nside,
        precompute_output_map=False,
        tensor_to_scalar=1e-3,
        map_dist=None,
        coord="C",
    ):
        """Websky CMB tensor-mode BB component

        Websky-compatible unlensed BB component due to primordial tensor perturbations
        The inputs are simulated with tensor-to-scalar ratio `r` of 1,
        then scaled by the `tensor_to_scalar` input parameter.

        Parameters
        ----------
        websky_version : str
            Websky version, see the documentation for more information
        nside : int
            Desired output HEALPix N_side
        precompute_output_map : bool
            If True, the output map is precomputed in the constructor
        tensor_to_scalar : float
            Tensor to scalar ratio `r`, ratio between the tensor and the
            scalar perturbations power spectra
        map_dist : pysm.MapDist
            see the PySM documentation
        """

        filename = utils.RemoteData(coord).get(
            "websky/{}/tensor_BB_r_1_cl.fits".format(websky_version)
        )
        self.tensor_to_scalar = tensor_to_scalar
        super().__init__(
            filename,
            input_units="uK_CMB",
            input_reference_frequency=None,
            nside=nside,
            precompute_output_map=precompute_output_map,
            has_polarization=True,
            from_cl=True,
            from_cl_seed=0,  # always do same realization
            map_dist=map_dist,
        )

    def compute_output_map(self, alm):
        return super().compute_output_map(alm) * np.sqrt(self.tensor_to_scalar)


class WebSkyCMB(PrecomputedAlms):
    def __init__(
        self,
        websky_version,
        nside,
        precompute_output_map=False,
        seed=1,
        lensed=True,
        map_dist=None,
        coord="C",
    ):
        filename = utils.RemoteData(coord).get(
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


class WebSkyCMBMap(pysm.CMBMap):
    def __init__(
        self,
        websky_version,
        nside,
        precompute_output_map=False,
        seed=1,
        lensed=True,
        include_solar_dipole=False,
        map_dist=None,
        coord="C",
    ):
        template_nside = 512 if nside <= 512 else 4096
        lens = "" if lensed else "un"
        soldip = "solardipole_" if include_solar_dipole else ""
        filenames = [
            utils.RemoteData(coord).get(
                f"websky/{websky_version}/map_{pol}_{lens}lensed_alm_seed{seed}_{soldip}nside{template_nside}.fits"
            )
            for pol in "IQU"
        ]
        super().__init__(
            map_I=filenames[0],
            map_Q=filenames[1],
            map_U=filenames[2],
            nside=nside,
            map_dist=map_dist,
        )
