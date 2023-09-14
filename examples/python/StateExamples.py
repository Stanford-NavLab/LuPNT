import pylupnt as pnt
import numpy as np

import sys
sys.path.append("pylupntutil")
import pntautodiff as ad

def degToRad(deg):
    return deg * np.pi / 180

x = ad.realvec([6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341])
cart = pnt.CartesianState(x, pnt.CoordSystem.ITRF)
print(" ")
print("Cartesian state:")
cart.print()

p = 11067.790
e = 0.83285
a = p / (1 - pow(e, 2))
i = degToRad(87.87)
Omega = degToRad(227.89)
w = degToRad(53.38)
nu = degToRad(92.335)
M = pnt.true_to_mean(nu, e)

x_oe = pnt.ClassicalOE(a, e, i, Omega, w, M, cs=pnt.CoordSystem.ITRF)
print(" ")
print("Classical orbital elements:")
x_oe.print()

print("")
print("a = ", x_oe.a().val())

# Convert orbital elements to Cartesian state
x_cart_from_oe = pnt.coe_to_cart(x_oe, mu=pnt.MU_EARTH)
print(" ")
print("Converted Cartesian State:")
x_cart_from_oe.print()
