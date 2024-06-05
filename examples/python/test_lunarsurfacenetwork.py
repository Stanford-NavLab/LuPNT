import numpy as np
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import pylupnt as pnt

import pylupnt.plots as plots
import pylupnt.utils as u

from pylupnt import SpiceInterface as sp


def degToRad(deg):
    return deg * pnt.RAD_PER_DEG

# check if folder exists, if not, create it
parent_dir = '/mnt/g/My Drive/Stanford_Documents/Documents/Research/Lunar_Surface_Station_ION2024'
sub_dir = 'simplebatch_keplarian'
show_plots = False

# Simulation parameters
rng_obs_std_dev_m = 0.1 #1   # (in meters) standard deviation of passive ranging measurements
init_pos_std = 0.1 #100 #km
init_vel_std = init_pos_std*1e-3

# dynamics models
truth_dynamics = pnt.KeplerianDynamics(pnt.MU_MOON)
batchfilter_dynamics = pnt.CartesianTwoBodyDynamics(pnt.MU_MOON, integrator='RK4')
# batchfilter_dynamics = pnt.NBodyDynamics()
# batchfilter_dynamics.set_primary_body(pnt.Body.MOON)
# batchfilter_dynamics.add_body(pnt.Body.EARTH)


# 
dt_meas = 60 # seconds
num_sat_periods = 5

# batch filter params
tol_batch = 1e-4
verbose_in_batch_filter = False
max_batch = 10


# Define satellite parameters (ELFO, Lunar Pathfinder)
# 13.2 hr ELFO, defined in a J2000 frame
a = 6539.1 # km
e = 0.60
i = degToRad(74.54)
Omega = degToRad(4.20)
w = degToRad(92.12)
M = degToRad(180)
tai_string = "2022/08/01 00:00:00"

# An ~11 hr ELFO, but defined in OP frame -- need to create conversion function
# a = 5760 # km
# e = 0.58
# i = degToRad(54.856)
# Omega = 0
# w = degToRad(86.322)
# M = degToRad(180)
# tai_string = "2030/10/01 00:00:00"




colors = plots.COLORS
fig_vis = plots.Plot3D(elev=-25, azim=-50, figsize=(8, 8))
fig_vis.plot_surface(pnt.MOON)
fig_vis.set_labels("x", "y", "z")
plt.title("Satellite orbit (PA)")


###########################################################################
#
# Create folder to save plots
#

# concatenate parent and sub directory
fldr_path = os.path.join(parent_dir, sub_dir)
# check if folder exists, if not, create it
run_num = 1
run_num_str = 'run' + str(run_num)
save_path = os.path.join(fldr_path, run_num_str)
while os.path.exists(save_path):
    print('run_num', run_num, 'already exists!')
    run_num += 1
    run_num_str = 'run' + str(run_num)
    save_path = os.path.join(fldr_path, run_num_str)
# create directory
os.makedirs(save_path)
print('Directory to save: ', save_path)
#
#
###########################################################################

def get_range(pos_arr1, pos_arr2):
    # assert that both arrays are of dimension X x 3
    assert pos_arr1.shape[1] == 3 and pos_arr2.shape[1] == 3

    return np.linalg.norm(pos_arr1 - pos_arr2, axis=1)

def get_range_elevation_and_visibility(pos_arr1, ref_pos_arr2, el_mask_deg=5):
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
    el_rad_arr = np.zeros((len(rng_arr)))
    for i in range(len(pos_arr1)):
        [az, el, rng] = pnt.cartesian_to_azimuth_elevation_range(ref_pos_arr2[i,:], pos_arr1[i,:])
        el_rad_arr[i] = el
        vis_arr[i] = (el >= el_mask_rad)
    return rng_arr, el_rad_arr, vis_arr

# Solve WLS problem (W is diagonal -- so only give diagonal elements as vector)
# def solve_wls(H, b, W=None):
#     # WLS solution
#     # H: measurement matrix of dim (n x s) -- n is num meas, s is state dim
#     # b: ranging observations of dim (n)
#     # W: weight vector of dim (n) (diagonal of weight matrix) -- if None, set to ones
#     # x = (H^T W H)^{-1} H^T W b
#     if W is None:
#         W = np.ones(H.shape[0])
    
#     # get H tilde and b tilde
#     sqrt_W = np.sqrt(W)
#     H_tilde = H * sqrt_W[:, np.newaxis]
#     b_tilde = b * sqrt_W

#     # solve for x
#     x = np.linalg.lstsq(H_tilde, b_tilde, rcond=None)[0]
#     return x

# get H matrix row 
def get_H_row(x_sat, x_lss):
    # x_cart: satellite position in inertial frame
    # x_stat: station position in inertial frame
    # H = [u/r, v/r, w/r, u_dot, v_dot, w_dot]
    r = np.linalg.norm(x_sat - x_lss)
    u = x_sat[0] - x_lss[0]
    v = x_sat[1] - x_lss[1]
    w = x_sat[2] - x_lss[2]
    return np.array([u/r, v/r, w/r, 0, 0, 0])


# specify initial time epoch (TAI frame)
tai_time = sp.string_to_tai(tai_string)
print('TAI time:', tai_time, type(tai_time))

oe = np.array([a, e, i, Omega, w, M])

# State
x_oe = pnt.ClassicalOE(oe, coord_sys=pnt.CoordSystem.MI) # MI -- J200 frame
print(" ")
print("Classical orbital elements:")
print(x_oe.vector)

# get orbital period from orbital elements
n = np.sqrt(pnt.MU_MOON / a ** 3)
T = 2*np.pi / n
print('Orbital period [hrs]:', T/3600)

# Add base station
az_stat_rad = pnt.RAD_PER_DEG * 0
el_stat_rad = pnt.RAD_PER_DEG * (-90) 
r_stat_km = pnt.R_MOON 
az_el_rkm_stat = np.array([az_stat_rad, el_stat_rad, r_stat_km])

# convert spherical coordinates to cartesian
theta = pnt.RAD_PER_DEG * 180
phi = 0
lss_pos_PA_km = pnt.spherical_to_cartesian(r_sph=[r_stat_km, theta, phi]).reshape([-1,3])
# concatenate r_stat_cart_km with zeros to get 1x6 vector
lss_pv_PA_km = np.concatenate((lss_pos_PA_km, np.zeros([1,3])), axis=1).reshape([-1,6])
print('Lunar station 6D state (Cartesian):', lss_pv_PA_km)
print()

# Propagate satellite orbit for certain amount of time and save state
tot_duration = num_sat_periods * T
t_arr = np.arange(0, tot_duration + dt_meas, dt_meas)
t_arr_tai = tai_time + t_arr 
len_t_arr = len(t_arr)
x_oe_arr = np.zeros((len(t_arr), 6))
sat_true_pv_MI_km_arr = np.zeros((len(t_arr), 6))
sat_true_pv_PA_km_arr = np.zeros((len(t_arr), 6))
lss_pv_MI_km_arr = np.zeros((len(t_arr), 6))
for i, t in enumerate(t_arr):
    truth_dynamics.propagate(x_oe, dt_meas)
    # save propagated state
    x_oe_arr[i,:] = x_oe.vector

    # convert from orbital elements to cartesian (inertial frame)
    # x_cart_arr[i,:] = pnt.classical_to_cartesian(x_oe, pnt.MU_MOON).vector
    sat_MI_pos = pnt.classical_to_cartesian(x_oe, pnt.MU_MOON).vector

    sat_true_pv_MI_km_arr[i,:] = sat_MI_pos 
    
    # convert base station position from PA frame to MI (at current TAI time)
    curr_tai_time = t_arr_tai[i]
    rv_stat_MI_km = pnt.CoordConverter.convert(curr_tai_time, lss_pv_PA_km, coord_sys_in=pnt.CoordSystem.PA, coord_sys_out=pnt.CoordSystem.MI)
    lss_pv_MI_km_arr[i,:] = rv_stat_MI_km

    # get satellite position in PA frame (just for plotting)
    sat_true_pv_PA_km_arr[i,:] = pnt.CoordConverter.convert(curr_tai_time, sat_MI_pos, coord_sys_in=pnt.CoordSystem.MI, coord_sys_out=pnt.CoordSystem.PA)


# plot satellite position (in PA frame)
fig_vis.plot(sat_true_pv_PA_km_arr[0:int(len_t_arr),0:3])

# Plot base station (in PA frame)
fig_vis.plot(lss_pos_PA_km, '*')
if show_plots:
    plt.show()

# 4. Compute the range and range rate from the satellite to the base station
rng_km_arr, el_rad_arr, vis_arr = get_range_elevation_and_visibility(sat_true_pv_MI_km_arr[:, 0:3], lss_pv_MI_km_arr[:, 0:3])
el_deg_arr = el_rad_arr / pnt.RAD_PER_DEG


# induce noise in the range measurements
Ri_km = np.array([[(rng_obs_std_dev_m * 1e-3)**2]])
inv_Ri_km = np.linalg.inv(Ri_km)
rng_obs_km_arr = rng_km_arr + np.random.normal(0, rng_obs_std_dev_m*1e-3, len(rng_km_arr))

# get times, ranges, and elevations when visible
t_arr_vis = t_arr[vis_arr == 1]
rng_arr_vis = rng_km_arr[vis_arr == 1]
el_deg_arr_vis = el_deg_arr[vis_arr == 1]

# plot ranges when visible
plt.figure()
plt.plot(t_arr_vis/3600.0, rng_arr_vis, '.')
plt.grid()
plt.xlabel('Time [hrs]')
plt.ylabel('Range [km]')
plt.title('Range to lunar station, when visible')
if show_plots:
    plt.show()

# plot elevation wrt station when visible
plt.figure()
plt.plot(t_arr_vis/3600.0, el_deg_arr_vis, '.')
plt.grid()
plt.xlabel('Time [hrs]')
plt.ylabel('Elevation [deg]')
plt.title('Elevation to lunar station')
if show_plots:
    plt.show()





##########################################################################################################################
#
# Perform batch filtering for orbit determination
#
##########################################################################################################################

# Assume some a priori estimate of initial satellite state
# create noisy vector of length 6, with first 3 having std dev of 1e-3 and last 3 having std dev of 1e-6
P0 = np.diag([init_pos_std**2, init_pos_std**2, init_pos_std**2, init_vel_std**2, init_vel_std**2, init_vel_std**2])
init_noise = np.concatenate((np.random.normal(0, init_pos_std, 3), np.random.normal(0, init_vel_std, 3)), axis=0)
init_sat_est_pv_MI_km = sat_true_pv_MI_km_arr[0,:] + init_noise
print('initial error magnitude [km]:', np.linalg.norm(init_noise[0:3]), np.linalg.norm(init_noise[3:6]))


# start with state deviation of 0
state_dev = np.zeros((6,1))
X0_star = init_sat_est_pv_MI_km.reshape([-1,1]) + state_dev

# things to plot later
norm_state_dev_change_pos_arr = []
norm_state_dev_change_vel_arr = []
init_sat_est_pv_MI_km_arr = np.zeros((len_t_arr, 6))
ave_pos_err_arr = []
ave_vel_err_arr = []

not_converged_yet = True
i_batch = 0
while i_batch < max_batch and not_converged_yet:

    # create array of cartesian estimates
    sat_est_pv_MI_km_arr = np.zeros((len_t_arr, 6))
    sat_est_pv_MI_km_arr[0,:] = X0_star.reshape([-1])

    # key initializations for batch filter
    Lambda = np.linalg.inv(P0)
    Eta = np.linalg.inv(P0)@state_dev    # used with Gamma to solve for updated state deviation
    Phi_0toim1 = np.eye(6)               # state transition matrix

    # start going forward through orbit, collecting observations, updating Gamma/Eta
    # iterator is i-1 (i minus 1, im1 for short), since we propagate to the next time step during the loop
    for im1 in range(len_t_arr-1):

        i = im1 + 1

        # using two-body dynamics, get the STM to convert to the next time step
        t_tai_im1 = t_arr_tai[im1]
        t_tai_i = t_arr_tai[i]
        dt_integ = (t_tai_i - t_tai_im1)/10.0

        # Integrate reference trajectory and STM from t_(i-1) to t_i
        # Propagation occurs in km
        sat_pv_im1 = sat_est_pv_MI_km_arr[i-1,:]
        if verbose_in_batch_filter and i==1:
            print()
            print('i =', i)
            print('State at i-1:', sat_pv_im1)
            print()
        sat_pv_i, Phi_im1toi = batchfilter_dynamics.propagate_with_stm(sat_pv_im1.reshape([-1,1]), t_tai_im1, t_tai_i, dt_integ)
        Phi_0toi = Phi_im1toi @ Phi_0toim1    # Get complete STM from 0 to i
        if verbose_in_batch_filter and i==1:
            print()
            print('i =', i)
            print('TAI times:', t_tai_im1, t_tai_i)
            print('State at i-1:', sat_pv_im1)
            print('State at i :,', sat_pv_i)
            print('diff norm: ', sat_pv_i - sat_pv_im1)
            print()
        
        # Check if visible at time i -- if so, accumulate current observations    
        if vis_arr[i] == 1:

            # get observation and expected range
            rng_obs_km = rng_obs_km_arr[i].reshape([-1,1])
            rng_exp_km = get_range(sat_pv_i[0:3].reshape([-1,3]), lss_pv_MI_km_arr[i,0:3].reshape([-1,3])).reshape([-1,1])

            # get H matrix row wrt the current state and convert to that of initial state using the STM
            Hi_tilde = get_H_row(sat_pv_i[0:3], lss_pv_MI_km_arr[i,0:3]).reshape([1,-1])
            yi = rng_obs_km - rng_exp_km
            Hi = np.dot(Hi_tilde, Phi_0toi)


            if verbose_in_batch_filter and i == 1 or (i == 10000) or i == len_t_arr :
                print()
                print('--------------------------------------------')
                print('i =', i)
                print('TAI times:', t_tai_im1, t_tai_i)
                print('State at i-1:', sat_pv_im1)
                print('State at i :,', sat_pv_i)
                print()
                print('Hi tilde', Hi_tilde)
                print('STM i-1 to i', Phi_im1toi)
                print('STM 0 to i', Phi_0toi)
                print('Hi', Hi)
                print('yi:', yi)
                print('curr state', sat_pv_i)
                print()
                print('shape Hi_tilde', Hi_tilde.shape)
                print('shape Hi', Hi.shape)
                print('shape of yi', rng_obs_km.shape, rng_exp_km.shape)
                print('--------------------------------------------')
                print()

            # update Gamma and Eta
            Lambda = Lambda + Hi.T @ (inv_Ri_km @ Hi)
            Eta = Eta + Hi.T @ (inv_Ri_km @ yi)
        
        # updates for next iteration
        sat_est_pv_MI_km_arr[i,:] = sat_pv_i # save propagated state
        Phi_0toim1 = Phi_0toi # save STM for next iteration

        # print progress of i periodically
        if i % 3000 == 0:
            print('Current iteration:', i, ' (out of', str(len_t_arr) +')')
        
    
    if verbose_in_batch_filter:
        # plot propagated satellite state (in MI)
        fig_vis = plots.Plot3D(elev=-25, azim=-50, figsize=(10, 10))
        fig_vis.plot_surface(pnt.MOON)
        fig_vis.set_labels("x", "y", "z")
        plt.title("Propagated satellite orbit (MI)")
        fig_vis.plot(sat_est_pv_MI_km_arr[0:int(len_t_arr),0:3])
        if show_plots:
            plt.show()

        # plot difference between true and propagated state for x, y, and z on a 3x1 subplot
        plt.figure()
        plt.subplot(3,1,1)
        plt.plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,0])
        plt.grid()
        plt.subplot(3,1,2)
        plt.plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,1])
        plt.grid()
        plt.subplot(3,1,3)
        plt.plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,2])
        plt.grid()
        plt.xlabel('Time [hrs]')
        plt.title('Propagated satellite orbit (MI)')
        if show_plots:
            plt.show()

        # plot difference between true and propagated state for x, y, and z on a 3x1 subplot
        plt.figure()
        plt.xlabel('Time [hrs]')
        plt.subplot(3,1,1)
        plt.plot(t_arr/3600.0, sat_true_pv_MI_km_arr[:,0] - sat_est_pv_MI_km_arr[:,0])
        plt.grid()
        plt.subplot(3,1,2)
        plt.plot(t_arr/3600.0, sat_true_pv_MI_km_arr[:,1] - sat_est_pv_MI_km_arr[:,1])
        plt.grid()
        plt.subplot(3,1,3)
        plt.plot(t_arr/3600.0, sat_true_pv_MI_km_arr[:,2] - sat_est_pv_MI_km_arr[:,2])
        plt.grid()
        if show_plots:
            plt.show()


    
    
    # Now that we've exited the loop, solve for the updated state deviation
    state_dev_hat = np.linalg.inv(Lambda) @ Eta

    # update the initial state estimate and stat deviation
    X0_star = X0_star + state_dev_hat
    state_dev = state_dev - state_dev_hat

    # update batch filter iteration and check convergence
    i_batch += 1
    norm_state_dev = np.linalg.norm(state_dev_hat)

    # append to position and velocity norm error for this iteration
    ave_pos_err_arr.append(np.linalg.norm(sat_true_pv_MI_km_arr[:,0:3] - sat_est_pv_MI_km_arr[:,0:3]))
    ave_vel_err_arr.append(np.linalg.norm(sat_true_pv_MI_km_arr[:,3:6] - sat_est_pv_MI_km_arr[:,3:6]))
    norm_state_dev_change_pos_arr.append(np.linalg.norm(state_dev_hat[0:3]))
    norm_state_dev_change_vel_arr.append(np.linalg.norm(state_dev_hat[3:6]))
    if norm_state_dev < tol_batch:
        print()
        print('Converged to solution? Yes!')
        print('     State deviation:', state_dev.T)
        print('     Final error vector:', X0_star.T-sat_true_pv_MI_km_arr[0,:])
        not_converged_yet = False
    else:
        print()
        print('Converged to solution? No')
        print('     State deviation:', state_dev)

    # if iteration is the first, update the initial array
    if i_batch == 0:
        # update with a deepy copy
        init_sat_est_pv_MI_km_arr = np.copy(sat_est_pv_MI_km_arr)

        
    
# Compare with true orbit
# from the updated X0_star, propagate the orbit and compare with true orbit
# create array of cartesian estimates
sat_est_pv_MI_km_arr = np.zeros((len_t_arr, 6))
sat_est_pv_MI_km_arr[0,:] = X0_star.reshape([-1])
for i in range(1, len_t_arr):
    # using two-body dynamics, get the STM to convert to the next time step
    t_tai_im1 = t_arr_tai[i-1]
    t_tai_i = t_arr_tai[i]
    dt_integ = (t_tai_i - t_tai_im1)/10.0

    # Integrate reference trajectory and STM from t_(i-1) to t_i
    sat_pv_im1 = sat_est_pv_MI_km_arr[i-1,:]
    sat_pv_i, Phi_im1toi = batchfilter_dynamics.propagate_with_stm(sat_pv_im1.reshape([-1,1]), t_tai_im1, t_tai_i, dt_integ)
    sat_est_pv_MI_km_arr[i,:] = sat_pv_i # save propagated state

# plot difference between true and propagated state for x, y, and z on a 3x1 subplot
# fig, axs = plt.subplots(3) 
# axs[0].plot(t_arr/3600.0, sat_true_pv_MI_km_arr[:,0] - sat_est_pv_MI_km_arr[:,0])
# axs[1].plot(t_arr/3600.0, sat_true_pv_MI_km_arr[:,1] - sat_est_pv_MI_km_arr[:,1])
# axs[2].plot(t_arr/3600.0, sat_true_pv_MI_km_arr[:,2] - sat_est_pv_MI_km_arr[:,2])
# fig.show()

fig, ax = plt.subplots(3,1)
ax[0].set_title('Difference between true and filter''s estimated orbit (MI)')
# .title('Difference between true and filter''s estimated orbit (MI)')
ax[0].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,0] - sat_true_pv_MI_km_arr[:,0])
ax[0].set_ylabel('x [km]')
ax[0].grid()
ax[1].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,1] - sat_true_pv_MI_km_arr[:,1])
ax[1].set_ylabel('y [km]')
ax[1].grid()
ax[2].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,2] - sat_true_pv_MI_km_arr[:,2])
ax[2].set_ylabel('z [km]')
ax[2].grid()
ax[2].set_xlabel('Time [hrs]')
plt.savefig(save_path + '/final_orbital_error_xyz.svg')
if show_plots:
    plt.show()

plt.figure()
plt.title('Total orbital error')
plt.plot(t_arr/3600.0, np.linalg.norm(sat_est_pv_MI_km_arr - sat_true_pv_MI_km_arr, axis=1))
plt.grid()
plt.xlabel('Time [hrs]')
plt.ylabel('Error [km]')
plt.savefig(save_path + '/final_orbital_error_total.svg')
if show_plots:
    plt.show()


fig, ax = plt.subplots(2,1)
ax[0].set_title('Iterative updates in state deviation vector')
ax[0].plot(np.arange(len(norm_state_dev_change_pos_arr))+1, norm_state_dev_change_pos_arr)
ax[0].set_ylabel('Norm position update [km]')
ax[0].grid()
ax[1].plot(np.arange(len(norm_state_dev_change_vel_arr))+1, norm_state_dev_change_vel_arr)
ax[1].set_ylabel('Norm velocity update [km/s]')
ax[1].grid()
ax[1].set_xlabel('batch filter iterations')
plt.savefig(save_path + '/update_over_iterations.svg')
if show_plots:
    plt.show()


fig, ax = plt.subplots(2,1)
ax[0].set_title('Average error in position and velocity')
ax[0].plot(np.arange(len(ave_pos_err_arr)), ave_pos_err_arr)
ax[0].set_ylabel('Position error [km]')
ax[0].grid()
ax[1].plot(np.arange(len(ave_vel_err_arr)), ave_vel_err_arr)
ax[1].set_ylabel('Velocity error [km/s]')
ax[1].grid()
ax[1].set_xlabel('batch filter iterations')
plt.savefig(save_path + '/error_over_iterations.svg')
if show_plots:
    plt.show()


# plot initial satellite orbit error
fig, ax = plt.subplots(3,1)
ax[0].set_title('Initial satellite position error')
ax[0].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,0] - sat_true_pv_MI_km_arr[:,0])
ax[0].set_ylabel('x [km]')
ax[0].grid()
ax[1].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,1] - sat_true_pv_MI_km_arr[:,1])
ax[1].set_ylabel('y [km]')
ax[1].grid()
ax[2].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,2] - sat_true_pv_MI_km_arr[:,2])
ax[2].set_ylabel('z [km]')
ax[2].grid()
ax[2].set_xlabel('Time [hrs]')
plt.savefig(save_path + '/init_orbit_error_pos.svg')
if show_plots:
    plt.show()


fig, ax = plt.subplots(3,1)
ax[0].set_title('Initial satellite velocity error')
ax[0].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,3] - sat_true_pv_MI_km_arr[:,3])
ax[0].set_ylabel('x vel [km/s]')
ax[0].grid()
ax[1].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,4] - sat_true_pv_MI_km_arr[:,4])
ax[1].set_ylabel('y vel [km/s]')
ax[1].grid()
ax[2].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,5] - sat_true_pv_MI_km_arr[:,5])
ax[2].set_ylabel('z vel [km/s]')
ax[2].grid()
ax[2].set_xlabel('Time [hrs]')
plt.savefig(save_path + '/init_orbit_error_vel.svg')
if show_plots:
    plt.show()


# if similar order of magnitude, plot initial/final on the same plot
fig, ax = plt.subplots(3,1)
ax[0].set_title('Initial/final satellite position error')
ax[0].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,0] - sat_true_pv_MI_km_arr[:,0])
ax[0].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,0] - sat_true_pv_MI_km_arr[:,0])
ax[0].legend(['Initial', 'Final'])
ax[0].set_ylabel('x [km]')
ax[0].grid()
ax[1].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,1] - sat_true_pv_MI_km_arr[:,1])
ax[1].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,1] - sat_true_pv_MI_km_arr[:,1])
ax[1].set_ylabel('y [km]')
ax[1].legend(['Initial', 'Final'])
ax[1].grid()
ax[2].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,2] - sat_true_pv_MI_km_arr[:,2])
ax[2].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,2] - sat_true_pv_MI_km_arr[:,2])
ax[2].set_ylabel('z [km]')
ax[2].legend(['Initial', 'Final'])
ax[2].grid()
ax[2].set_xlabel('Time [hrs]')
plt.savefig(save_path + '/init_final_orbit_error_pos.svg')
if show_plots:
    plt.show()



fig, ax = plt.subplots(3,1)
ax[0].set_title('Initial/final satellite velocity error')
ax[0].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,3] - sat_true_pv_MI_km_arr[:,3])
ax[0].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,3] - sat_true_pv_MI_km_arr[:,3])
ax[0].set_ylabel('x vel [km/s]')
ax[0].legend(['Initial', 'Final'])
ax[0].grid()
ax[1].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,4] - sat_true_pv_MI_km_arr[:,4])
ax[1].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,4] - sat_true_pv_MI_km_arr[:,4])
ax[1].set_ylabel('y vel [km/s]')
ax[1].legend(['Initial', 'Final'])
ax[1].grid()
ax[2].plot(t_arr/3600.0, init_sat_est_pv_MI_km_arr[:,5] - sat_true_pv_MI_km_arr[:,5])
ax[2].plot(t_arr/3600.0, sat_est_pv_MI_km_arr[:,5] - sat_true_pv_MI_km_arr[:,5])
ax[2].set_ylabel('z vel [km/s]')
ax[2].legend(['Initial', 'Final'])
ax[2].grid()
ax[2].set_xlabel('Time [hrs]')
plt.savefig(save_path + '/init_final_orbit_error_vel.svg')
if show_plots:
    plt.show()



