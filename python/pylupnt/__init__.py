try:
    from .pylupnt_pybind import *  # py2 py3 compatible

    from .utils import *
    from . import plot
    from . import render
    from .math_utils import *
    from . import scenarios
    from . import crater_detection

except ImportError as e:
    # this was installed with as a python wheel
    from pylupnt_pybind import *

    import plot
    import render
    from utils import *
    from math_utils import *
    import scenarios
    import dataset

try:
    import pkg_resources  # part of setuptools

    __version__ = pkg_resources.require("pylupnt")[0].version
except Exception:
    __version__ = "@PROJECT_VERSION@"
