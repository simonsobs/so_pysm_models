import numpy as np

from .. import get_so_models
import pytest

components_dict = {comp[0]: comp for comp in ["dust", "synchrotron", "freefree", "ame"]}

# Expected I and Q emission in uK_RJ at pixel 100
expected = {
    "SO_a1": [0.0004612855020969098, 4.015659454536007e-06],
    "SO_s1": [0.3903231775760541, -0.024661686876035677],
    "SO_d1": [1.5400970036635495, 0.044328112549527524],
    "SO_a0": [0.0004949370177550887, 0.0],
    "SO_f0": [0.2165784407714049, 0.0],
    "SO_s0": [0.4133453059990284, -0.02504213606690909],
    "SO_d0": [1.9386241336315098, 0.05124213250215371],
}


@pytest.mark.parametrize(
    "model_tag", ["SO_d0", "SO_s0", "SO_f0", "SO_a0", "SO_d1", "SO_s1", "SO_a1"]
)
def test_get_so_models(model_tag):
    """Test the `get_so_models` function

    Only check that we can access the 512 version of the template
    and that the result has no NaN"""
    from pysm import Sky

    component_name = components_dict[model_tag.split("_")[1][0]]
    sky = Sky({component_name: get_so_models(model_tag, nside=128)})
    component = getattr(sky, component_name)
    emission = component(nu=100)
    assert not np.any(np.isnan(emission))
    # Compare I and Q at pixel 100
    for IQ in [0, 1]:
        assert emission[IQ][100] == pytest.approx(expected[model_tag][IQ])
