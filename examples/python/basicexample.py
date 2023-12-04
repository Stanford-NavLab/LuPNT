import numpy as np
from pylupnt import pylupnt_pybind as pnt

from inspect import getmembers, isfunction, isclass, ismodule

print("Versions:", pnt.__version__)
print("File:", pnt.__file__)

print("List of Classes")
print(getmembers(pnt, isclass))
print(" ")

print("List of Functions")
print(getmembers(pnt, isfunction))

print("List of Modules")
print(getmembers(pnt, ismodule))