import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylupnt as pnt

import pylupnt.plots as plots
import pylupnt.utils as u

# mpl.rcParams.update(mpl.rcParamsDefault)

colors = plots.COLORS
# 3d plot
fig_vis = plots.Plot3D(elev=-25, azim=-50, figsize=(10, 10))
fig_vis.plot_surface(pnt.MOON)
fig_vis.set_labels("x", "y", "z")
plt.title("Satellite orbit (MI)")
# print(colors)

def degToRad(deg):
    print(pnt.RAD_PER_DEG)
    return deg * pnt.RAD_PER_DEG

def get_range(pos_arr1, pos_arr2):
    # assert that both arrays are of dimension X x 3
    assert pos_arr1.shape[1] == 3 and pos_arr2.shape[1] == 3

    diff_vect = pos_arr1 - pos_arr2
    return np.linalg.norm(pos_arr1 - pos_arr2, axis=1)

def get_range_and_visibility(pos_arr1, ref_pos_arr2, el_mask_deg=5):
    # convert elevation mask to radians
    el_mask_rad = pnt.RAD_PER_DEG * el_mask_deg

    # assert that both arrays are of dimension X x 3
    assert pos_arr1.shape[1] == 3 and ref_pos_arr2.shape[1] == 3

    # compute the range vector
    diff_vect = pos_arr1 - ref_pos_arr2
    # compute the range
    rng_arr = np.linalg.norm(diff_vect, axis=1)
    # # compute the unit vector
    # unit_vect = diff_vect / rng[:, np.newaxis]

    # convert unit vector to elevation, azimuth, range
    vis_arr = np.zeros((len(rng_arr)))
    for i in range(len(pos_arr1)):
        [az, el, rng] = pnt.cartesian_to_azimuth_elevation_range(ref_pos_arr2, pos_arr1[i,:])
        vis_arr[i] = (el >= el_mask_rad)
    return rng_arr, vis_arr


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
    # TODO: coordinate transformation here? from inertial to moon-centered inertial


fig_vis.plot(x_cart_arr[0:int(len_t_arr/2),0:3])

# Add base station
az_stat_rad = pnt.RAD_PER_DEG * 0
el_stat_rad = pnt.RAD_PER_DEG * (-90) 
r_stat_km = pnt.R_MOON 
az_el_rkm_stat = np.array([az_stat_rad, el_stat_rad, r_stat_km])

# convert spherical coordinates to cartesian
theta = pnt.RAD_PER_DEG * 180
phi = 0
r_stat_cart_km = pnt.spherical_to_cartesian(r_sph=[r_stat_km, theta, phi]).reshape([-1,3])
print('Lunar station at:', r_stat_cart_km)
print()
fig_vis.plot(r_stat_cart_km, '*')
plt.show()

# 4. Compute the range and range rate from the satellite to the base station
rng_arr, vis_arr = get_range_and_visibility(x_cart_arr[:, 0:3], r_stat_cart_km)

print(rng_arr.shape)
print(t_arr.shape)

# Use matplotlib to plot the range and use latex for the labels
plt.figure()
plt.plot(t_arr, rng_arr)
plt.grid()
plt.xlabel('Time [s]')
plt.ylabel('Range [km]')
plt.title('Range to lunar station')
plt.show()

plt.figure()
plt.plot(t_arr, vis_arr)
plt.grid()
plt.xlabel('Time [s]')
plt.ylabel('Visibility')
plt.title('Visibility of lunar station')
plt.show()

# get times and ranges when visible
t_arr_vis = t_arr[vis_arr == 1]
rng_arr_vis = rng_arr[vis_arr == 1]
# print(len(t_arr_vis))
# print('sum of vis arr', np.sum(vis_arr))

# plot ranges when visible
plt.plot()
plt.plot(t_arr_vis, rng_arr_vis)
plt.grid()
plt.xlabel('Time [s]')
plt.ylabel('Range [km]')
plt.title('Range to lunar station, when visible')
plt.show()


