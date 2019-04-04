from . import InterpolatingComponent
from . import utils
import pysm

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

    def get_filenames(self,path):
        """Get filenames for a websky version

        For a standard interpolating component, we list files in folder,
        here we need to know the names in advance so that we can only download the required maps
        """
        websky_version = path
        if websky_version == "0.3":
            filenames = {}
            def add_freq(frequency_GHz): filenames[frequency_GHz] = 'websky/0.3/cib_'+str(frequency_GHz).zfill(4)+ '.fits'
            for base_freq in [27, 39, 93, 145, 225, 280]:
                for delta_freq in [-1, 0, 1]:
                    add_freq(base_freq + delta_freq)
            for base_freq in [100,217,353,545,857]:
                delta_freq = 0
                add_freq(base_freq + delta_freq)
        return filenames


    def read_map(self, freq):
        filename = utils.get_data_from_url(self.maps[freq])
        return self.read_map_file(freq, filename)
