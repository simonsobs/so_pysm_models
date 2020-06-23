# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This sub-module is destined for common non-package specific utility
# functions.

import os
from astropy.utils import data
import warnings

DATAURL = {
    "C": "https://portal.nersc.gov/project/cmb/so_pysm_models_data/equatorial/",
    "G": "https://portal.nersc.gov/project/cmb/so_pysm_models_data/",
}
PREDEFINED_DATA_FOLDERS = {
    "C": [
        "/global/project/projectdirs/cmb/www/so_pysm_models_data/equatorial",  # NERSC
        "/mnt/sdceph/users/zonca/simonsobs/so_pysm_models_data/equatorial",  # SDSC
    ],
    "G": [
        "/global/project/projectdirs/cmb/www/so_pysm_models_data/",  # NERSC
        "/mnt/sdceph/users/zonca/simonsobs/so_pysm_models_data/",  # SDSC
    ],
}


class RemoteData:
    def __init__(self, coord):
        """Access template from remote server

        PySM input templates are stored on the CMB project space at NERSC
        and are made available via web.
        The get method of this class tries to access data locally from one
        of the PREDEFINED_DATA_FOLDERS defined above, if it fails, it
        retrieves the files and caches them remotely using facilities
        provided by `astropy.utils.data`.

        Parameters
        ----------

        coord : "C" or "G"
            We have templates either in Equatorial or Galactic reference frame
        """
        self.data_url = DATAURL[coord]
        self.data_folders = PREDEFINED_DATA_FOLDERS[coord]

    def get(self, filename):
        for folder in self.data_folders:
            full_path = os.path.join(folder, filename)
            if os.path.exists(full_path):
                warnings.warn(f"Access data from {full_path}")
                return full_path
        with data.conf.set_temp("dataurl", self.data_url), data.conf.set_temp(
            "remote_timeout", 30
        ):
            warnings.warn(f"Retrieve data for {filename} (if not cached already)")
            map_out = data.get_pkg_data_filename(filename, show_progress=True)
        return map_out
