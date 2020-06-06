import numpy as np
import healpy as hp
from astropy.utils.data import get_pkg_data_filename

from .. import GaussianDust
from astropy.tests.helper import assert_quantity_allclose

try:  # PySM >= 3.2.1
    import pysm3.units as u
except ImportError:
    import pysm.units as u


def test_gaussian_dust():
    test_map_filename = get_pkg_data_filename(
        "data/Gaussian_dust_353GHz_uRJ_testmap_ns128.fits.zip"
    )
    test_map = hp.read_map(test_map_filename, (0, 1, 2)) << u.uK_RJ

    nside = hp.npix2nside(len(test_map[0]))
    gaussian_dust = GaussianDust(nside=nside, seed=18)
    m = gaussian_dust.get_emission(353*u.GHz)

    assert_quantity_allclose(m, test_map)
