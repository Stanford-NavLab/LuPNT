import numpy as np

# Classical Orbital Elements
# (a, e, i, Omega, w, M) [m, -, rad]
# - Elliptical Lunar Frozen Orbit
coe_array_elfo = np.array(
    [9750.5, 0.7, np.deg2rad(63.5), np.deg2rad(90.0), np.deg2rad(0.0), np.deg2rad(30.0)]
)
# - Low Lunar Orbit
coe_array_llo = np.array(
    [100.0, 1e-3, np.deg2rad(28.5), np.deg2rad(0.0), np.deg2rad(0.0), np.deg2rad(20.0)]
)

# Quasi-Nonsingular Relative Orbital Elements
# (ada, adl, adex, adey, adix, adiy) [m]
roe_array_1 = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
roe_array_2 = np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0])

# Print options
array2string_kwargs = {
    "separator": ", ",
    "formatter": {"float_kind": lambda x: "%g" % x},
}


def unpack_gmat(rvector6):
    return np.array([rvector6.Get(i) for i in range(6)])
