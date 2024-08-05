"""Lunar Positioning, Navigation, and Timing (LuPNT) Library"""

# here we import the contents of our compiled C++ module
try:
    from .pylupnt import *  # py2 py3 compatible
    from .python_code import pure_python_list
except Exception:
    # this was installed with as a python wheel
    from pylupnt import *
    from python_code import pure_python_list
