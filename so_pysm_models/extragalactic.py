from . import InterpolatingComponent


class WebSkyCIB(InterpolatingComponent):
    """PySM component interpolating between precomputed maps"""

    def __init__(
        self,
        path,
        input_units="MJysr",
        target_nside=4096,
        interpolation_kind="linear",
        pixel_indices=None,
        mpi_comm=None,
        verbose=False,
    ):
        super().__init__(
            self,
            path,
            input_units,
            target_nside,
            interpolation_kind,
            False,
            pixel_indices,
            mpi_comm,
            verbose,
        )

    def get_fnames(self, path):
        # Override this to implement name convention
        from websky_model import websky

        ws = websky.WebSky(path)
        freqs = [100, 143, 217, 353, 545]
        fnames = {}
        for freq in freqs:
            fnames[freq] = ws.cib_map_file_name(str(freq), "0")
        return fnames
