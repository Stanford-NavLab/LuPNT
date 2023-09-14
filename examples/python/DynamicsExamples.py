import pylupnt as pnt
import numpy as np

import sys

sys.path.append("pylupntutil")
import pntautodiff as ad


def degToRad(deg):
    return deg * np.pi / 180


p = 11067.790
e = 0.83285
a = p / (1 - pow(e, 2))
i = degToRad(87.87)
Omega = degToRad(227.89)
w = degToRad(53.38)
nu = degToRad(92.335)
M = pnt.true_to_mean(nu, e)

# State
x_oe = pnt.ClassicalOE(a, e, i, Omega, w, M, cs=pnt.CoordSystem.GCRF)
print(" ")
print("Classical orbital elements:")
x_oe.print()

# Dynamics
mu_moon = 4902.800066
keplerian_dynamics = pnt.KeplerianDynamics(mu_moon)
print(" ")
print("Keplerian dynamics object:")
print(keplerian_dynamics)

# Propagate
t_end = 10.0
keplerian_dynamics.propagate(x_oe, t_end)
print(" ")
print("Propagated state:")
x_oe.print()
