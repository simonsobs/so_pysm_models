import os.path

import numpy as np
import healpy as hp

import pytest

import pysm

from .. import PrecomputedAlms


@pytest.fixture()
def setup(tmpdir):
    # tmpdir is a py.test feature to provide a temporary folder

    folder = tmpdir.mkdir("alms")

    np.random.seed(12)
    alm_size = hp.Alm.getsize(lmax=100)
    alms = 1j * np.random.normal(size=(3, alm_size))
    alms += np.random.normal(size=(3, alm_size))

    # str needed to support Python 3.5
    filename = os.path.join(str(folder), "alms.fits")
    hp.write_alm(filename, alms)
    return alms, filename


def test_precomputed_alms(setup):

    alms, filename = setup

    nside = 64
    # we assume the original `alms` are in `K_CMB`
    ref_freq = 40
    test_map_K_CMB = hp.alm2map(alms, nside=nside)

    alms_K_RJ = alms * pysm.convert_units("K_CMB", "K_RJ", ref_freq)
    filename_K_RJ = filename.replace(".fits", "_RJ.fits")
    hp.write_alm(filename_K_RJ, alms_K_RJ)

    precomputed_alms = PrecomputedAlms(
        filename=filename_K_RJ,
        nside=nside,
        input_units="uK_RJ",
        input_reference_frequency_GHz=ref_freq,
    )
    m = precomputed_alms.signal(23)

    np.testing.assert_allclose(
        m, test_map_K_CMB * pysm.convert_units("uK_CMB", "uK_RJ", 23)
    )

    freqs = np.array([1, 10, 100])
    m_multifreq = precomputed_alms.signal(freqs)

    assert m_multifreq.shape == (3, 3, hp.nside2npix(64))

    for freq, m in zip(freqs, m_multifreq):
        np.testing.assert_allclose(
            m, test_map_K_CMB * pysm.convert_units("uK_CMB", "uK_RJ", freq)
        )


def test_precomputed_alms_K_CMB(setup):

    alms, filename = setup

    nside = 64
    test_map = hp.alm2map(alms, nside=nside)
    precomputed_alms = PrecomputedAlms(
        filename=filename, nside=nside, input_units="K_CMB"
    )

    freqs = np.array([1, 10, 100])
    m_multifreq = precomputed_alms.signal(freqs)

    assert m_multifreq.shape == (3, 3, hp.nside2npix(64))

    for freq, m in zip(freqs, m_multifreq):
        np.testing.assert_allclose(
            m, test_map * pysm.convert_units("K_CMB", "uK_RJ", freq)
        )
