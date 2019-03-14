# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This sub-module is destined for common non-package specific utility
# functions.

import os
from astropy.utils import data

DATAURL = "http://portal.nersc.gov/project/cmb/so_pysm_models_data/"
PREDEFINED_DATA_FOLDERS = ["/global/project/projectdirs/cmb/www/so_pysm_models_data"]


def get_data_from_url(filename):
    """Retrieves input templates from remote server,
    in case data is available in one of the PREDEFINED_DATA_FOLDERS defined above,
    e.g. at NERSC, those are directly returned."""
    for folder in PREDEFINED_DATA_FOLDERS:
        full_path = os.path.join(folder, filename)
        if os.path.exists(full_path):
            return full_path
    with data.conf.set_temp("dataurl", DATAURL), data.conf.set_temp(
        "remote_timeout", 30
    ):
        map_out = data.get_pkg_data_filename(filename, show_progress=True)
    return map_out
