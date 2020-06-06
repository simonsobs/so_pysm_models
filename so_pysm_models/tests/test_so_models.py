import numpy as np

try:  # PySM >= 3.2.1
    import pysm3 as pysm
    import pysm3.units as u
except ImportError:
    import pysm.units as u
    import pysm

from .. import get_so_models
from astropy.tests.helper import assert_quantity_allclose
import pytest

components_dict = {comp[0]: comp for comp in ["dust", "synchrotron", "freefree", "ame"]}

# Expected I and Q emission in uK_RJ at pixel 98969
expected = {
    "SO_a1": [2.7826194464298837, 0.025349571461346275],
    "SO_s1": [11.109877, -0.020188138587129197],
    "SO_d1": [789.259705, 3.7430164427988215],
    "SO_a0": [1.4868411735277658, 0.0],
    "SO_f0": [4525.092680419471, 0.0],
    "SO_s0": [11.756732686207846, -0.021363560253508863],
    "SO_d0": [921.4452792523148, 4.369898779380189],
}


@pytest.mark.parametrize(
    "model_tag", ["SO_d0", "SO_s0", "SO_f0", "SO_a0", "SO_d1", "SO_s1", "SO_a1"]
)
def test_get_so_models(model_tag):
    """Test the `get_so_models` function

    Only check that we can access the 512 version of the template
    and that the result has no NaN"""

    sky = pysm.Sky(
        nside=128, component_objects=[get_so_models(model_tag, nside=128, coord="G")]
    )
    emission = sky.get_emission(freq=100 * u.GHz)

    assert not np.any(np.isnan(emission))
    # Compare I and Q at pixel 100
    for IQ in [0, 1]:
        if expected[model_tag][IQ] != 0:
            assert_quantity_allclose(
                emission[IQ][98969], expected[model_tag][IQ] * u.uK_RJ, rtol=1e-4
            )
