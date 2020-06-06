import healpy as hp
from astropy.utils.data import get_pkg_data_filename
try:  # PySM >= 3.2.1
    import pysm3.units as u
except ImportError:
    import pysm.units as u

from .. import GaussianSynchrotron
from astropy.tests.helper import assert_quantity_allclose


def test_gaussian_synchrotron():
    test_map_filename = get_pkg_data_filename(
        "data/Gaussian_synch_23GHz_uRJ_testmap_ns128.fits.zip"
    )
    test_map = hp.read_map(test_map_filename, (0, 1, 2)) << u.uK_RJ

    nside = hp.npix2nside(len(test_map[0]))
    gaussian_synch = GaussianSynchrotron(nside=nside, seed=31)
    m = gaussian_synch.get_emission(23*u.GHz)

    assert_quantity_allclose(m, test_map)
