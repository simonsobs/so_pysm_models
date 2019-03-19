from . import InterpolatingComponent
import os

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
            path,
            input_units,
            target_nside,
            interpolation_kind,
            has_polarization=False,
            pixel_indices=pixel_indices,
            mpi_comm=mpi_comm,
            verbose=verbose,
        )

    def get_filenames(self, path):
        filenames = {}
        for f in os.listdir(path):
            if f.endswith(".fits"):
                freq = int(f.split(".")[0].split("nu")[1])
                filenames[freq] = os.path.join(path, f)
        return filenames
