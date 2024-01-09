try:
    from .pylupnt_pybind import *  # py2 py3 compatible

    from . import plots, utils
except ImportError:
    # this was installed with as a python wheel
    from pylupnt_pybind import *

    from . import plots, utils

try:
    import pkg_resources  # part of setuptools

    __version__ = pkg_resources.require("pylupnt")[0].version
except Exception:
    __version__ = "@PROJECT_VERSION@"
