import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylupnt as pnt

import pylupnt.plots as plots
import pylupnt.utils as u

# mpl.rcParams.update(mpl.rcParamsDefault)

colors = plots.COLORS
# 3d plot
fig = plots.Plot3D(elev=-25, azim=-50, figsize=(10, 10))
fig.plot_surface(pnt.MOON)
fig.set_labels("x", "y", "z")
plt.title("Satellite orbit (MI)")
# print(colors)

def degToRad(deg):
    print(pnt.RAD_PER_DEG)
    return deg * pnt.RAD_PER_DEG

# Define satellite parameters (ELFO, Lunar Pathfinder)
a = 5760 # km
e = 0.58
i = degToRad(54.856)
Omega = 0
w = degToRad(86.322)
M = degToRad(180)
# nu = degToRad(92.335)
# M = pnt.true_to_mean_anomaly(nu, e)

oe = np.array([a, e, i, Omega, w, M])

# State
x_oe = pnt.ClassicalOE(oe, coord_sys=pnt.CoordSystem.MI)
print(" ")
print("Classical orbital elements:")
print(x_oe.vector)

# 2. Keplarian dynamics
keplarian_dynamics = pnt.KeplerianDynamics(pnt.MU_MOON)
print(keplarian_dynamics)

# get orbital period from orbital elements
n = np.sqrt(pnt.MU_MOON / a ** 3)
T = 2*np.pi / n
print('Orbital period [hrs]:', T/3600)

print("MU Earth", pnt.MU_EARTH)

# Propagate for certain amount of time and save propagated satellite state
dt = 1
t_arr = np.arange(0, 2*T + dt, dt)
len_t_arr = len(t_arr)
x_oe_arr = np.zeros((len(t_arr), 6))
x_cart_arr = np.zeros((len(t_arr), 6))
for i, t in enumerate(t_arr):
    keplarian_dynamics.propagate(x_oe, dt)
    # save propagated state
    x_oe_arr[i,:] = x_oe.vector

    # convert from orbital elements to cartesian
    x_cart_arr[i,:] = pnt.classical_to_cartesian(x_oe, pnt.MU_MOON).vector

fig.plot(x_cart_arr[0:int(len_t_arr/2),0:3])
# plt.show()

# Add base station
az_stat_rad = pnt.RAD_PER_DEG * 0
el_stat_rad = pnt.RAD_PER_DEG * (-90) 
r_stat_km = pnt.R_MOON 
az_el_rkm_stat = np.array([az_stat_rad, el_stat_rad, r_stat_km]);

# convert spherical coordinates to cartesian
theta = pnt.RAD_PER_DEG * 180
phi = 0
r_stat_cart_km = pnt.spherical_to_cartesian(r_sph=[r_stat_km, theta, phi])
print('Lunar station at:', r_stat_cart_km)
print()
fig.plot(r_stat_cart_km.reshape([-1,3]), '*')
plt.show()

# 4. Compute the range and range rate from the satellite to the base station

