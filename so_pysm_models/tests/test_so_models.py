import numpy as np

from .. import get_so_models
import pytest

components_dict = {comp[0]: comp for comp in ["dust", "synchrotron", "freefree", "ame"]}


@pytest.mark.parametrize("model_tag", ["SO_d0", "SO_s0", "SO_f0", "SO_a0"])
def test_get_so_models(model_tag):
    """Test the `get_so_models` function

    Only check that we can access the 512 version of the template
    and that the result has no NaN"""
    from pysm import Sky

    component_name = components_dict[model_tag.split("_")[1][0]]
    sky = Sky({component_name: get_so_models(model_tag, nside=128)})
    component = getattr(sky, component_name)
    assert not np.any(np.isnan(component(nu=100)))
