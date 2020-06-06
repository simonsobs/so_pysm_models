import numpy as np
import healpy as hp

from astropy.utils.data import get_pkg_data_filename
try:  # PySM >= 3.2.1
    import pysm3.units as u
except ImportError:
    import pysm.units as u

from .. import PrecomputedAlms
from astropy.tests.helper import assert_quantity_allclose


def test_precomputed_alms():
    alms_filename = get_pkg_data_filename(
        "data/fullskyUnlensedUnabberatedCMB_alm_set00_00000.fits.zip"
    )
    save_name = get_pkg_data_filename("data/test_cmb_map.fits.zip")
    nside = 32
    # Make an IQU sim
    precomputed_alms = PrecomputedAlms(
        alms_filename, nside=nside, input_units="uK_CMB",
        has_polarization=True,
        #input_reference_frequency=148*u.GHz
    )
    simulated_map = precomputed_alms.get_emission(148*u.GHz).to(u.uK_CMB, equivalencies=u.cmb_equivalencies(148*u.GHz))
    expected_map = hp.read_map(save_name, field=(0, 1, 2)) << u.uK_CMB
    assert_quantity_allclose(simulated_map, expected_map)
    assert simulated_map.shape[0] == 3
