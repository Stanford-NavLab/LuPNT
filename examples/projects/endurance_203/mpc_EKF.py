import os
import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import pylupnt as pnt
import scipy as sc
import cvxpy as cvx
from time import time
from tqdm.auto import tqdm

np.random.seed(11)

class EKF:
    def __init__(self, x0, P0, lat_long_0, elev_grid, heading, cell_elevation) -> None:
        # define the process noise covariance
        # self.Q = np.eye(3) * 0.0001
        # define the measurement noise covariance
        # self.R = np.eye(3) * 0.0001
        # define the initial state estimate
        self.x_hat = x0
        # define the initial state estimate covariance
        self.P = P0
        # self.s_mpc = 
        self.lat_long_0 = lat_long_0
        self.elev_grid = elev_grid
        self.cell_elevation = cell_elevation
        self.heading = heading
        self.res = 0.1


    def predict(self, control):
        # predict the next state estimate
        # print(f'Current EKF state: {self.x_hat}')
        # print(f'Cell elevation: {self.cell_elevation}')
        # print(f'Heading: {self.heading}')

        mpc_state = MoonPA_to_DEM(self.x_hat, self.cell_elevation, self.lat_long_0, self.heading)
        # print(f'Current MPC state: {mpc_state}')
        A, B, _ = affinize(mpc_state, control)
        mpc_state_next = dynamics(mpc_state, control)
        # print(f'Predicted state: {mpc_state_next}')
        self.x_hat, self.cell_elevation, self.heading = state2MoonPA(mpc_state_next, self.res, self.lat_long_0, self.elev_grid)
        Q =  np.eye(3) * 0.0001
        self.P = A @ self.P @ A.T + Q
        # return x_hat, P
    
    def update(self, satpos, true_state):
        # update the state estimate
        # true state and measurement from satellites
        n_meas = satpos.shape[0]
        los = np.zeros((n_meas, 3))
        for i in range(n_meas):
            los[i] = true_state - satpos[i]
        geom_range = np.linalg.norm(los, axis=1) # 1xN
        range_std = (5)/1000              # 5m    
        # range_std = 0
        measurement = geom_range + np.random.normal(0.0, range_std, size=n_meas)

        # measurement model and estimated parameters
        los_est = np.zeros((n_meas, 3))
        for i in range(n_meas):
            los_est[i] = self.x_hat - satpos[i]
        g_x = np.linalg.norm(los_est, axis=1)
        H = los_est / np.tile(g_x, (3, 1)).T
        R = np.eye(n_meas)* (range_std**2)

        # Kalman gain
        K = self.P @ H.T @ np.linalg.inv(H @ self.P @ H.T + R)
        #TODO: need measurement models
        self.x_hat = self.x_hat + K @ (measurement - g_x)
        self.P = (np.eye(3) - K @ H) @ self.P

    def ekf_step(self, control, satpos_t, true_state):

        # given a state and control, predict the next state
        # ekf.ekf_step(u_mpc[t, 0], satpos_t, s_true)

        # predicts in the MoonPA frame of reference
        self.predict(control)

        # true_state needs to be in PA frame
        self.update(satpos_t, true_state)

        # return self.x_hat
        




def dynamics(s, u, dt = 0.1):
    x, y, theta = s
    v, omega = u
    x_next = x + (dt * v * np.cos(theta))
    y_next = y + (dt * v * np.sin(theta))
    theta_next = theta + (dt * omega)
    s_next = np.array([x_next, y_next, theta_next])
    return s_next

def affinize(s, u, dt = 0.1):
    x, y, theta = s
    v, omega = u
    A = np.array([[1, 0, -dt * v * np.sin(theta)],
                        [0, 1, dt * v * np.cos(theta)],
                        [0, 0, 1]])
    B = dt * np.array([[np.cos(theta), 0],
                        [np.sin(theta), 0],
                        [0, 1]]) 
    c = dynamics(s, u) - A @ s - B @ u
    return A, B, c

def scp_iteration(f, s0, s_goal, s_prev, u_prev, P, Q, R, dt = 0.1, v_bound=5.0, omega_bound=3.0):
    n = s_prev.shape[-1]  # state dimension
    m = u_prev.shape[-1]  # control dimension
    N = u_prev.shape[0]  # number of steps
    # v_bound = 5 # km/hr (upper bound on rover velocity)
    # omega_bound = 3.0 # rad/hr ?
    
    A = np.zeros((N, n, n))
    B = np.zeros((N, n, m))
    c = np.zeros((N, n))
    
    # print(f' Previous traj: {s_prev}')
    
    for k in range(N):
        A_k, B_k, c_k = affinize(s_prev[k], u_prev[k], dt)
        A[k,:,:] = A_k
        B[k,:,:] = B_k
        c[k,:] = c_k
    
    s_cvx = cvx.Variable((N + 1, n))
    u_cvx = cvx.Variable((N, m))
    
    objective = 0.0
    constraints = []
    for k in range(N):
        objective += cvx.quad_form((s_cvx[k] - s_goal), Q) + cvx.quad_form(u_cvx[k], R) # sum the cost
        constraints.append(s_cvx[k+1] == A[k] @ s_cvx[k] + B[k] @ u_cvx[k] + c[k]) # dynamics constraint
        constraints.append(u_cvx[k, 0] <= v_bound)
        constraints.append(u_cvx[k, 0] >= 0.0)
        constraints.append(cvx.abs(u_cvx[k, 1]) <= omega_bound)
    objective += cvx.quad_form((s_cvx[N] - s_goal), P) # sum the terminal cost
    constraints.append(s_cvx[0] == s0) # add initial constraint
    # print(f's0 = {s0}')
    # print(f's_goal = {s_goal}')
    prob = cvx.Problem(cvx.Minimize(objective), constraints)
    prob.solve()
    if prob.status != "optimal":
        raise RuntimeError("SCP solve failed. Problem status: " + prob.status)
    s = s_cvx.value
    u = u_cvx.value
    J = prob.objective.value
    # print(f'u = {u[0]}')
    # print(f's_cvx[0] = {s[0]}')
    # print(f's_cvx[1] (linear) = {A[0] @ s[0] + B[0] @ u[0]}')
    # print(f's_cvx[1] = {s[1]}')
    # print(f'A = \n{A[0]}')
    # print(f'B = \n{B[0]}')
    return s, u, J

def solve_mpc(
    f,
    s0,
    s_goal,
    N,
    P,
    Q,
    R,
    eps,
    max_iters,
    dt=0.1,
    v_bound=5.0,
    omega_bound=3.0,
    s_init=None,
    u_init=None,
    convergence_error=False,
):
    """Solve the obstacle avoidance problem via SCP."""
    n = Q.shape[0]  # state dimension
    m = R.shape[0]  # control dimension


    s = np.linspace(s0, s_goal, N+1)
    u = np.ones((N, m))

    # print(f'Intial traj :{s}')
    
    # Do SCP until convergence or maximum number of iterations is reached
    converged = False
    J = np.zeros(max_iters + 1)
    J[0] = np.inf
    for i in range(max_iters):
        s, u, J[i + 1] = scp_iteration(f, s0, s_goal, s, u, P, Q, R, dt, v_bound, omega_bound)
        dJ = np.abs(J[i + 1] - J[i])
        if dJ < eps:
            converged = True
            # print(f'Converged in {i} iterations!')
            break
    if not converged and convergence_error:
        raise RuntimeError("SCP did not converge!")
    return s, u

def run_ekf_mpc(s0, s_goal, N, P, Q, R, T, N_scp, dt, v_bound, omega_bound, sat_data, tspan, t0, lat_long_0, elev_grid, sigma_init=10/1000):

    # sat_data (N_satellite, time, 6)

    n = 3 # state dimension
    m = 2 # control dimension
    eps = 1e-3 # SCP convergence tolerance

    f = dynamics
    s_mpc = np.zeros((T, N + 1, n))
    u_mpc = np.zeros((T, N, m))
    s_true = np.copy(s0)
    total_time = time()
    total_control_cost = 0.0
    s_init = None
    u_init = None

    # initialize the EKF as well (with some random noise added to the initial state estimate)

    # this is in the MoonPA frame
    # sigma_init = 10/1000 # 10m
    x0_hat, cell_elevation, heading = state2MoonPA(s0, 0.1, lat_long_0, elev_grid)
    # print(f'Check: {MoonPA_to_DEM(x0_hat, cell_elevation, lat_long_0, heading)}')
    x0_hat += np.random.randn(3) * sigma_init
    P0 = np.eye(3) * (sigma_init**2)
    ekf = EKF(x0_hat, P0, lat_long_0, elev_grid, heading, cell_elevation)
    t_true = t0
    s_est = MoonPA_to_DEM(x0_hat, ekf.cell_elevation, ekf.lat_long_0, ekf.heading)

    print(f'Initial state estimate: {ekf.x_hat}')
    print(f'Initial rover estimate: {s_est}')
    print(f'Intial covariance: {ekf.P}')
    print(f'Starting time: {t_true}')
    # print(f'Current EKF state: {ekf.x_hat}')
    # print(f'Cell elevation: {ekf.cell_elevation}')
    # print(f'Heading: {ekf.heading}')
    # print('----------')

    true_state_array = np.zeros((T, 3))
    est_state_array = np.zeros((T, 3))
    # true_state_array[0] = s_true


    for t in tqdm(range(T)):
        # solving MPC with the true starting state
        # s_mpc[t], u_mpc[t] = solve_mpc(f, s_true, s_goal, N, P, Q, R, eps, N_scp, dt, v_bound, omega_bound, s_init, u_init)

        # solving MPC with the state from the EKF
        s_mpc[t], u_mpc[t] = solve_mpc(f, s_est, s_goal, N, P, Q, R, eps, N_scp, dt, v_bound, omega_bound, s_init, u_init)

        # the next true state
        s_true = f(s_true, u_mpc[t, 0])
        x_true, cell_elevation, heading = state2MoonPA(s_true, 0.1, lat_long_0, elev_grid)
        # get the true time to match with simulation
        t_true += (3600*dt)
        t_idx = np.where(tspan == t_true)[0]

        # get the satellite data at this time
        satpos_t = sat_data[:, t_idx, :3]

        # the estimated state
        ekf.ekf_step(u_mpc[t, 0], satpos_t, x_true)
        s_est = MoonPA_to_DEM(ekf.x_hat, ekf.cell_elevation, ekf.lat_long_0, ekf.heading)

        # print(s_est)/
        # print(f'Estimated state: {s_est}')
        # print(f'True state: {s_true}')

        total_control_cost += u_mpc[t, 0].T @ R @ u_mpc[t, 0]

        # Use this solution to warm-start the next iteration
        u_init = np.concatenate([u_mpc[t, 1:], u_mpc[t, -1:]])
        s_init = np.concatenate(
            [s_mpc[t, 1:], f(s_mpc[t, -1], u_mpc[t, -1]).reshape([1, -1])]
        )

        true_state_array[t] = s_true
        est_state_array[t] = s_est

    total_time = time() - total_time
    # print("Total elapsed time:", total_time, "seconds")
    # print("Total control cost:", total_control_cost)
    # print(s_mpc[:, 0, 1])
    print(f'End time: {t_true}')
    return s_mpc, u_mpc, t_true, true_state_array, est_state_array


def plot_mpc(s0, s_goal, s_mpc, u_mpc, N, T, N_scp, n_waypt, dt):
    fig, ax = plt.subplots(2, 2, dpi=150, figsize=(8, 6))
    # fig.suptitle("$N = {}$, ".format(N) + r"$N_\mathrm{SCP} = " + "{}$".format(N_scp))
    fig.suptitle('State and Control Input Over Time')

    ms = 2
    # for t in range(T):
        # ax[0, 0].plot(s_mpc[t, :, 1], s_mpc[t, :, 0], "--*", color="k", markersize=ms, linewidth=0.5)
    ax[0, 0].plot(s_mpc[:, 0, 1], s_mpc[:, 0, 0], "-o", markersize=ms, linewidth=0.5)
    ax[0, 0].set_xlabel(r"$y(t)$")
    ax[0, 0].set_ylabel(r"$x(t)$")
    ax[0, 0].axis("equal")
    ax[0, 0].scatter(s0[1], s0[0], color='r', zorder=3,s=ms)
    ax[0, 0].scatter(s_goal[1], s_goal[0], color='r', zorder=3, s=ms)

    t_vec = np.linspace(0, T*(n_waypt - 1)*dt, len(s_mpc[:, 0, 2]))
    ax[0, 1].plot(t_vec, s_mpc[:, 0, 2]*180/np.pi, "-o", markersize=ms) #, label=r"$\theta(t)$")
    ax[0, 1].set_xlabel(r"$t$ [hr]")
    ax[0, 1].set_ylabel(r"$\theta(t)$ [deg]")
    ax[0, 1].set_xlim([0, t_vec[-1]])
    # ax[0, 1].legend()
    # ax[0, 1].axhline(s_goal[2]*180/np.pi)
    
    ls = 1
    for i in range(n_waypt-1):
        index = T * i
        ax[0, 1].axvline(t_vec[index],linestyle='--',color='k',zorder=0, linewidth=ls)
        ax[1, 0].axvline(t_vec[index],linestyle='--',color='k',zorder=0, linewidth=ls)
        ax[1, 1].axvline(t_vec[index],linestyle='--',color='k',zorder=0, linewidth=ls)
    ax[0, 1].axvline(t_vec[-1],linestyle='--',color='k',zorder=0, linewidth=ls)
    ax[1, 0].axvline(t_vec[-1],linestyle='--',color='k',zorder=0, linewidth=ls)
    ax[1, 1].axvline(t_vec[-1],linestyle='--',color='k',zorder=0, linewidth=ls)

    ax[1, 0].plot(t_vec, u_mpc[:, 0, 0], "-o", markersize=ms, linewidth=ls) #, label=r"$u_1(t)$")
    ax[1, 0].set_xlabel(r"$t$ [hr]")
    ax[1, 0].set_ylabel(r"$v(t)$ [km/hr]")
    ax[1, 0].set_xlim([0, t_vec[-1]])
    # ax[1, 0].legend()
    # ax[1, 0].set_ylim([-0.1, 5.1])

    ax[1, 1].plot(t_vec, u_mpc[:, 0, 1]*180/np.pi/3600, "-o", markersize=ms, linewidth=ls) #, label=r"$u_2(t)$")
    ax[1, 1].set_xlabel(r"$t$ [hr]")
    ax[1, 1].set_ylabel(r"$\omega(t)$ [deg/sec]")
    ax[1, 1].set_xlim([0, t_vec[-1]])
    # ax[1, 1].legend()
    # ax[1, 1].set_ylim(np.array([-3.1, 3.1])*180/np.pi)
    
    ax[0, 0].grid(True)

    # suffix = "_N={}_Nscp={}".format(N, N_scp)
    plt.tight_layout()
    # plt.savefig("soln_obstacle_avoidance" + suffix + ".png", bbox_inches="tight")
    plt.savefig("figures/mpc_state_control_over_time.png", bbox_inches="tight",dpi=300)

    # plt.show()



def state2MoonPA(state, res, lat_long_0, elev_grid):
    x, y, theta = state
    # print(f'x = {x}, y = {y}')
    x_grid = x/res; y_grid = y/res

    # print(f'x_grid = {x_grid}, y_grid = {y_grid}')

    # snap to the grid and get the corresponding indices
    state_idx = [int(x_grid), int(y_grid)]
    # get the latitude and longitude of the current state
    R_Moon = 1737.4  # km
    lat = lat_long_0[0] + x/R_Moon
    long = lat_long_0[1] + y/R_Moon
    cell_elevation = elev_grid[state_idx[0], state_idx[1]]

    R_Moon = 1737.4  # km
    user_mcmf_pos_x = (R_Moon + cell_elevation) * np.cos(lat) * np.cos(long)
    user_mcmf_pos_y = (R_Moon + cell_elevation) * np.cos(lat) * np.sin(long)
    user_mcmf_pos_z = (R_Moon + cell_elevation) * np.sin(lat)

    return np.array([user_mcmf_pos_x, user_mcmf_pos_y, user_mcmf_pos_z]), cell_elevation, theta

def latlong_to_MoonPA(lat, long, cell_elevation):
    # convert latitude and longitude to MoonPA
    # cell_elevation is a float
    R_Moon = 1737.4  # km
    user_mcmf_pos_x = (R_Moon + cell_elevation) * np.cos(lat) * np.cos(long)
    user_mcmf_pos_y = (R_Moon + cell_elevation) * np.cos(lat) * np.sin(long)
    user_mcmf_pos_z = (R_Moon + cell_elevation) * np.sin(lat)

    # 3x1 position vector
    user_mcmf_loc = np.array([user_mcmf_pos_x, user_mcmf_pos_y, user_mcmf_pos_z])
    return user_mcmf_loc

def MoonPA_to_DEM(xyz_MoonPA, elev, lat_long_0, heading):
    # given x, y, z on a sphere, find the latitude and longitude
    x, y, z = xyz_MoonPA
    R_Moon = 1737.4  # km
    lat = np.arcsin(z / (R_Moon+elev))
    long = np.arctan2(y, x)

    x_mpc = R_Moon*(lat - lat_long_0[0])
    y_mpc = R_Moon*(long - lat_long_0[1])

    return np.array([x_mpc, y_mpc, heading])

    


