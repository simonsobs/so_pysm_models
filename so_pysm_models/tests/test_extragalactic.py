import numpy as np
import healpy as hp

from .. import WebSkyCIB


def test_cib(tmp_path):

    nside = 4
    shape = hp.nside2npix(nside)
    hp.write_map(tmp_path / "nu0010.fits", np.zeros(shape, dtype=np.float32))
    hp.write_map(tmp_path / "nu0020.fits", np.ones(shape, dtype=np.float32))

    interp = WebSkyCIB(
        tmp_path, "uK_RJ", nside, interpolation_kind="linear",
    )

    interpolated_map = interp.signal(nu=15)
    np.testing.assert_allclose(0.5 * np.ones((1, shape)), interpolated_map)
