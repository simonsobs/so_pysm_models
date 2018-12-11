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

    filename = os.path.join(folder, "alms.fits")
    hp.write_alm(filename, alms)
    return alms, filename


def test_precomputed_alms(setup):

    alms, filename = setup

    nside = 64
    test_map = hp.alm2map(alms, nside=nside)
    precomputed_alms = PrecomputedAlms(
        filename=filename, target_nside=nside, input_units="uK_RJ"
    )
    m = precomputed_alms.signal(23)

    np.testing.assert_allclose(m, test_map)

    m_multifreq = precomputed_alms.signal(np.array([1, 10, 100]))

    assert m_multifreq.shape == (3, 3, hp.nside2npix(64))

    for each in m_multifreq:
        np.testing.assert_allclose(each, test_map)


def test_precomputed_alms_K_CMB(setup):

    alms, filename = setup

    nside = 64
    test_map = hp.alm2map(alms, nside=nside)
    precomputed_alms = PrecomputedAlms(
        filename=filename, target_nside=nside, input_units="K_CMB"
    )

    freqs = np.array([1, 10, 100])
    m_multifreq = precomputed_alms.signal(freqs)

    assert m_multifreq.shape == (3, 3, hp.nside2npix(64))

    for freq, m in zip(freqs, m_multifreq):
        np.testing.assert_allclose(
            m, test_map * pysm.convert_units("K_CMB", "uK_RJ", freq)
        )
