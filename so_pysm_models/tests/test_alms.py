import os.path

import numpy as np
import healpy as hp

from .. import PrecomputedAlms


def test_gaussian_synchrotron(tmpdir):
    # tmpdir is a py.test feature to provide a temporary folder

    folder = tmpdir.mkdir("alms")

    np.random.seed(12)
    alm_size = hp.Alm.getsize(lmax=100)
    alms = 1j * np.random.normal(size=(3, alm_size))
    alms += np.random.normal(size=(3, alm_size))

    filename = os.path.join(folder, "alms.fits")
    hp.write_alm(filename, alms)

    nside = 64

    test_map = hp.alm2map(alms, nside=nside)
    precomputed_alms = PrecomputedAlms(filename=filename, target_nside=nside)
    m = precomputed_alms.signal(23)

    np.testing.assert_allclose(m, test_map)
