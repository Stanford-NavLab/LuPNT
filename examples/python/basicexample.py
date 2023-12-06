import numpy as np
import pylupnt as pnt

from inspect import getmembers, isfunction, isclass, ismodule

# print("Versions:", pnt.__version__)
# print("File:", pnt.__file__)

print("List of Classes")
print(getmembers(pnt, isclass))
print(" ")

print("List of Modules")
print(getmembers(pnt, ismodule))

temp = pnt.wrapToPi(3.14)
print(temp)
