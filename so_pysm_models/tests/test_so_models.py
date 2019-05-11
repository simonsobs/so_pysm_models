import numpy as np

from .. import get_so_models
import pytest

components_dict = {comp[0]: comp for comp in ["dust", "synchrotron", "freefree", "ame"]}

# Expected I and Q emission in uK_RJ at pixel 98969
expected = {
    "SO_a1": [2.7826194464298837,0.025349571461346275],
    "SO_s1": [9.512234698510888,-0.020188138587129197],
    "SO_d1": [747.7160343440288,3.7430164427988215],
    "SO_a0": [1.4868411735277658,0.0],
    "SO_f0": [4525.092680419471,0.0],
    "SO_s0": [11.756732686207846,-0.021363560253508863],
    "SO_d0": [921.4452792523148,4.369898779380189],
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
        assert emission[IQ][98969] == pytest.approx(expected[model_tag][IQ])
