import pylupnt as pnt
import numpy as np


def deg2Rad(deg):
    return deg * np.pi / 180


p = 11067.790
e = 0.83285
a = p / (1 - pow(e, 2))
i = deg2Rad(87.87)
Omega = deg2Rad(227.89)
w = deg2Rad(53.38)
nu = deg2Rad(92.335)
M = pnt.true2mean_anomaly(nu, e)

oe = np.array([a, e, i, Omega, w, M])

# State
x_oe = pnt.ClassicalOE(oe, frame=pnt.Frame.MOON_CI)
print(" ")
print("Classical orbital elements:")
print(x_oe.vector)

# Dynamics
mu_moon = 4902.800066
keplerian_dynamics = pnt.KeplerianDynamics(mu_moon)

# Propagate
t_end = 10.0
keplerian_dynamics.propagate(x_oe, t_end)
print(" ")
print("Propagated state:")
print(x_oe.vector)
