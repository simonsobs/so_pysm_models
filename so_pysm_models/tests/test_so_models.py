import numpy as np

from .. import get_so_models


def test_get_so_models():
    """Test the `get_so_models` function

    Only check that we can access the 512 version of the template
    and that the result has no NaN"""
    from pysm import Sky
    sky = Sky({
            "dust" : get_so_models("SO_d0", nside=128),
    })
    assert not np.any(np.isnan(sky.dust(nu=100)))
