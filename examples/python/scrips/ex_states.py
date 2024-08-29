import pylupnt as pnt
import numpy as np


def deg2Rad(deg):
    return deg * np.pi / 180


x_cart = np.array([6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341])
print(" ")
print("Cartesian state:")
print(x_cart)

p = 11067.790
e = 0.83285
a = p / (1 - pow(e, 2))
i = deg2Rad(87.87)
Omega = deg2Rad(227.89)
w = deg2Rad(53.38)
nu = deg2Rad(92.335)
M = pnt.true2mean_anomaly(nu, e)

x_oe = pnt.ClassicalOE(np.array([a, e, i, Omega, w, M]), frame=pnt.Frame.ITRF)
print(" ")
print("Classical orbital elements:")
print(x_oe.vector)

print("")
print("a = ", x_oe.a)

# Convert orbital elements to Cartesian state
x_cart_from_oe = pnt.classical2cart(x_oe, mu=pnt.GM_EARTH)
print(" ")
print("Converted Cartesian State:")
print(x_cart_from_oe.vector)
