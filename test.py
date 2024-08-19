import pylupnt as pnt
import numpy as np

dyn = pnt.KeplerianDynamics(pnt.GM_EARTH)

# Classical orbital elements (a, e, i, W, w, M) [km, -, rad, rad, rad, rad]
coe0 = np.array([5740, 0.58, 54.856, 0, 86.322, 0])
coe0[2:] *= pnt.RAD

a = coe0[0]
t0 = pnt.gregorian2time(2022, 1, 15, 0, 0, 0)
T = 2 * np.pi * np.sqrt(np.power(a, 3) / pnt.GM_EARTH)
tspan = np.linspace(0, T, 20)
tfs = t0 + tspan

coes = dyn.propagate(coe0, t0, tfs[0])
print(coes)
# help(pnt.KeplerianDynamics)
