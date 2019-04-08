import numpy as np
import healpy as hp

from .. import utils
from .. import WebSkyCIB, WebSkySZ


def test_cib(tmp_path, monkeypatch):

    monkeypatch.setattr(utils, "PREDEFINED_DATA_FOLDERS", [str(tmp_path)])
    nside = 4
    shape = hp.nside2npix(nside)

    path = tmp_path / "websky" / "0.3"
    path.mkdir(parents=True)
    hp.write_map(path / "cib_0094.fits", np.zeros(shape, dtype=np.float32))
    hp.write_map(path / "cib_0100.fits", np.ones(shape, dtype=np.float32))

    interp = WebSkyCIB("0.3", "uK_RJ", nside, interpolation_kind="linear")

    interpolated_map = interp.signal(nu=97)
    np.testing.assert_allclose(
        np.interp(97, [94, 100], [0, 1]) * np.ones((1, shape)), interpolated_map
    )

def test_ksz(tmp_path, monkeypatch):

    monkeypatch.setattr(utils, "PREDEFINED_DATA_FOLDERS", [str(tmp_path)])
    nside = 4
    shape = hp.nside2npix(nside)

    path = tmp_path / "websky" / "0.3"
    path.mkdir(parents=True)
    hp.write_map(path / "ksz.fits", np.ones(shape, dtype=np.float32))
    hp.write_map(path / "cib_0100.fits", np.ones(shape, dtype=np.float32))

    ksz = WebSkySZ("0.3", sz_type="kinetic", target_nside=nside)

    ksz_map = ksz.signal(nu=100)
    np.testing.assert_allclose(
        np.ones(ksz_map.shape)*0.777228, ksz_map, rtol=1e-4
    )

def test_tsz(tmp_path, monkeypatch):

    monkeypatch.setattr(utils, "PREDEFINED_DATA_FOLDERS", [str(tmp_path)])
    nside = 4
    shape = hp.nside2npix(nside)

    path = tmp_path / "websky" / "0.3"
    path.mkdir(parents=True)
    hp.write_map(path / "tsz.fits", np.ones(shape, dtype=np.float32)*1e-6)
    hp.write_map(path / "cib_0100.fits", np.ones(shape, dtype=np.float32))

    tsz = WebSkySZ("0.3", sz_type="thermal", target_nside=nside)

    tsz_map = tsz.signal(nu=100)
    np.testing.assert_allclose(
        np.ones(tsz_map.shape)*-3.193671, tsz_map, rtol=1e-4
    )
