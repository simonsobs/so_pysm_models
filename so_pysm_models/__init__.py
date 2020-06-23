# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *

# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__minimum_python_version__ = "3.5"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple(
    (int(val) for val in __minimum_python_version__.split("."))
):
    raise UnsupportedPythonError(
        "so_pysm_models does not support Python < {}".format(__minimum_python_version__)
    )

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    pass

from .synchrotron import GaussianSynchrotron
from .dust import GaussianDust
from .alms import PrecomputedAlms
from .co_lines import COLines
from .extragalactic import WebSkySZ, WebSkyCIB, WebSkyCMB, WebSkyCMBMap, WebSkyCMBTensor
from .so_models import get_so_models
