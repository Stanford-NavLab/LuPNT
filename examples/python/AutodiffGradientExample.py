import pylupnt as pnt
import numpy as np

import sys
sys.path.append("pylupntutil")
import pntautodiff as ad

# function defintions
def myfunc(x):
    y = ad.dotvec(x, ad.expvec(x)) + ad.dotvec(ad.cosvec(x), ad.sinvec(x))
    return y

# main script
x_val = [1.0, 2.0, 3.0]
x = ad.realvec(x_val)
print("x = ")
print(x)
print(" ")

F = myfunc(x)
print("F = ")
print(F)
print(" ")

g = pnt.gradient(myfunc, x)
print("Gradient Vector = ")
print(g)
print(" ")

# numerical gradient
grad_num = ad.numerical_gradient(myfunc, x)

print("Numerical Gradient Vector =")
print(grad_num)
print(" ")