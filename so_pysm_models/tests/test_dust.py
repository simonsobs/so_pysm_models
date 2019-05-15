import numpy as np
import healpy as hp
from astropy.utils.data import get_pkg_data_filename

from .. import GaussianDust


def test_gaussian_dust():
    test_map_filename = get_pkg_data_filename(
        "data/Gaussian_dust_353GHz_uRJ_testmap_ns128.fits.zip"
    )
    test_map = hp.read_map(test_map_filename, (0, 1, 2))

    nside = hp.npix2nside(len(test_map[0]))
    gaussian_dust = GaussianDust(nside=nside, seed=18)
    m = gaussian_dust.signal(353)

    np.testing.assert_allclose(m, test_map)
