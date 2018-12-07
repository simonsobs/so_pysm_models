import numpy as np
import healpy as hp
from astropy.utils.data import get_pkg_data_filename

from .. import GaussianSynchrotron


def test_gaussian_synchrotron():
    test_map_filename = get_pkg_data_filename(
        "data/Gaussian_synch_23GHz_uRJ_testmap.fits.zip"
    )
    test_map = hp.read_map(test_map_filename, (0, 1, 2))

    nside = hp.npix2nside(len(test_map[0]))
    gaussian_synch = GaussianSynchrotron(target_nside=nside, seed=31)
    m = gaussian_synch.signal(23)

    np.testing.assert_allclose(m, test_map)
