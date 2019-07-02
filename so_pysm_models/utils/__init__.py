# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This sub-module is destined for common non-package specific utility
# functions.

import os
from astropy.utils import data
import warnings

DATAURL = "http://portal.nersc.gov/project/cmb/so_pysm_models_data/equatorial"
PREDEFINED_DATA_FOLDERS = [
    "/global/project/projectdirs/cmb/www/so_pysm_models_data/equatorial",  # NERSC
    "/simons/scratch/zonca/simonsobs/so_pysm_models_data/equatorial",  # SDSC
]


def get_data_from_url(filename):
    """Retrieves input templates from remote server,
    in case data is available in one of the PREDEFINED_DATA_FOLDERS defined above,
    e.g. at NERSC, those are directly returned."""
    for folder in PREDEFINED_DATA_FOLDERS:
        full_path = os.path.join(folder, filename)
        if os.path.exists(full_path):
            warnings.warn(f"Access data from {full_path}")
            return full_path
    with data.conf.set_temp("dataurl", DATAURL), data.conf.set_temp(
        "remote_timeout", 30
    ):
        warnings.warn(f"Retrieve data for {filename} (if not cached already)")
        map_out = data.get_pkg_data_filename(filename, show_progress=True)
    return map_out
