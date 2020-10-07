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
        "data/Planck_bestfit_alm_seed_583_lmax_95_K_CMB.fits.zip"
    )
    save_name = get_pkg_data_filename("data/map_nside_32_from_Planck_bestfit_alm_seed_583_K_CMB.fits.zip")
    nside = 32
    # Make an IQU sim
    precomputed_alms = PrecomputedAlms(
        alms_filename, nside=nside, input_units="uK_CMB",
        has_polarization=True,
        #input_reference_frequency=148*u.GHz
    )
    simulated_map = precomputed_alms.get_emission(148*u.GHz).to(u.uK_CMB, equivalencies=u.cmb_equivalencies(148*u.GHz))
    expected_map = hp.read_map(save_name, field=(0, 1, 2)) << u.uK_CMB
    assert simulated_map.shape[0] == 3
    assert_quantity_allclose(simulated_map, expected_map)

def test_precomputed_alms_clip():
    # If clipping to lmax of 95 were not applied, this test would fail
    alms_filename = get_pkg_data_filename(
        "data/Planck_bestfit_alm_seed_583_lmax_120_K_CMB.fits.zip"
    )
    save_name = get_pkg_data_filename("data/map_nside_32_from_Planck_bestfit_alm_seed_583_K_CMB.fits.zip")
    nside = 32
    # Make an IQU sim
    precomputed_alms = PrecomputedAlms(
        alms_filename, nside=nside, input_units="uK_CMB",
        has_polarization=True,
        #input_reference_frequency=148*u.GHz
    )
    simulated_map = precomputed_alms.get_emission(148*u.GHz).to(u.uK_CMB, equivalencies=u.cmb_equivalencies(148*u.GHz))
    expected_map = hp.read_map(save_name, field=(0, 1, 2)) << u.uK_CMB
    assert simulated_map.shape[0] == 3
    assert_quantity_allclose(simulated_map, expected_map)
