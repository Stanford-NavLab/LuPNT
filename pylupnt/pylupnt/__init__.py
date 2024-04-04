try:
    from .pylupnt_pybind import *  # py2 py3 compatible

    from .utils import *
    from .plots import *
except ImportError:
    # this was installed with as a python wheel
    from pylupnt_pybind import *

    from utils import *
    from plots import *

try:
    import pkg_resources  # part of setuptools

    __version__ = pkg_resources.require("pylupnt")[0].version
except Exception:
    __version__ = "@PROJECT_VERSION@"
