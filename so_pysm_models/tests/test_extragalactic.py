import numpy as np
import healpy as hp

from .. import utils
from ..import WebSkyCIB


def test_cib(tmp_path, monkeypatch):

    monkeypatch.setattr(utils, "PREDEFINED_DATA_FOLDERS",  [str(tmp_path)])
    nside = 4
    shape = hp.nside2npix(nside)

    path = tmp_path / "websky" / "0.2"/"cib"
    path.mkdir(parents=True)
    hp.write_map(path / "cib_ns4096_nu0094.fits", np.zeros(shape, dtype=np.float32))
    hp.write_map(path / "cib_ns4096_nu0144.fits", np.ones(shape, dtype=np.float32))

    interp = WebSkyCIB(
        "0.2", "uK_RJ", nside, interpolation_kind="linear",
    )

    interpolated_map = interp.signal(nu=100)
    np.testing.assert_allclose(np.interp(100, [94,144], [0,1]) * np.ones((1, shape)), interpolated_map)
