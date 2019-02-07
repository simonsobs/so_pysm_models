import numpy as np
import healpy as hp

from astropy.utils.data import get_pkg_data_filename

from .. import PrecomputedAlms


def test_precomputed_alms():
    alms_filename = get_pkg_data_filename(
        "data/fullskyUnlensedUnabberatedCMB_alm_set00_00000.fits.zip"
    )
    save_name = get_pkg_data_filename("data/test_cmb_map.fits.zip")
    nside = 32
    # Make an IQU sim
    precomputed_alms = PrecomputedAlms(
        alms_filename, target_nside=nside, input_units="uK_RJ",
        input_reference_frequency_GHz=148, has_polarization=True
    )
    simulated_map = precomputed_alms.signal()
    expected_map = hp.read_map(save_name, field=(0, 1, 2))
    np.testing.assert_allclose(simulated_map, expected_map)
    assert simulated_map.shape[0] == 3
