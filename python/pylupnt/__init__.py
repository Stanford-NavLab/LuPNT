"""An example of how to use pybind11 and numpy"""

# here we import the contents of our compiled C++ module
from ._pylupnt import *

# we can also import from python modules as usual:
from .python_code import pure_python_list
