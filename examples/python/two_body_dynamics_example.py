import pylupnt as pnt
import numpy as np


def degToRad(deg):
    return deg * np.pi / 180

p = 11067.790
e = 0.83285
a = p / (1 - pow(e, 2))
i = degToRad(87.87)
Omega = degToRad(227.89)
w = degToRad(53.38)
nu = degToRad(92.335)
M = pnt.true_to_mean_anomaly(nu, e)

oe = np.array([a, e, i, Omega, w, M])

# State
x_oe = pnt.ClassicalOE(oe, coord_sys=pnt.CoordSystem.MI)
# print(" ")
# print("Classical orbital elements:")
# print(x_oe.vector)

x_cart = pnt.classical_to_cartesian(x_oe, pnt.MU_MOON).vector

# Dynamics
mu_moon = 4902.800066
cart_dynamics = pnt.CartesianTwoBodyDynamics(mu_moon)

# Propagate
t_end = 600.0
dt = t_end / 10.0
print(' ')
print('Initial state:')
print(x_cart)
x_prop = cart_dynamics.propagate(x_cart, 0.0, t_end, dt)
print(" ")
print("Propagated state:")
print(x_prop)

