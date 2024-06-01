import os
import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import pylupnt as pnt
import scipy as sc
import cvxpy as cvx


class Rover_Agent:
    def __init__(self, rover_lat_long):
        # constants
        self.g = 1.62       # m/s^2
        self.m = 487        # kg (Endurance report states 487 fully margined)
        self.dt = 0.01      # s
        self.N = 50         # grid size
        self.R_moon = 1737.4e3  # m

        # rover location
        self.rover_loc = latlong_to_MoonPA(rover_lat_long[0], rover_lat_long[1])

        # define any relevant control constraints here
        # define any additional rover constants here


    def dynamics(self, state, control):
        """Differential Drive Robot (continuous inertial dynamics)"""
        # 3 state variables --> x, y, theta
        # 2 control inputs --> v, omega
        x, y, theta = state
        v, omega = control
        # return dynamics
        return jnp.array([v * jnp.cos(theta), v * jnp.sin(theta), omega])
    
    def ss_model_dynamics(self, state, control):
        """Discrete time state space model"""
        # Assumption --> Euler integration
        # x_t+1 = x_t + dt * f(x_t, u_t) where f is the dynamics function
        x, y, theta = state
        v, omega = control

        # define the state space model
        A = self.dt * jnp.array([[1, 0, -v * jnp.sin(theta)],
                        [0, 1, v * jnp.cos(theta)],
                        [0, 0, 1]])
        
        B = self.dt * jnp.array([[jnp.cos(theta), 0],
                        [jnp.sin(theta), 0],
                        [0, 1]]) 
        
        # return the state space model
        return A, B
    
    # def ss_model_measurement(self, state):
        # """Discrete time measurement model"""
        # # Assumption --> perfect measurements
        # # y_t = x_t + w_t where w_t is the measurement noise
        # # define the measurement model
        # H = jnp.eye(3)
        # # return the measurement model
        # return H

        # consider the case where the rover has an onboard receiver for nav satellites

    def propagate_dynamics(self, state, control):
        """Propagate the dynamics of the
        rover using the state space model"""
        A, B = self.ss_model_dynamics(state, control)
        # returns the next state
        return A @ state + B @ control
    
def latlong_to_MoonPA(lat, long):
    # convert latitude and longitude to MoonPA
    R_Moon = 1737.4  # m
    user_mcmf_pos_x = R_Moon * np.cos(lat) * np.cos(long)
    user_mcmf_pos_y = R_Moon * np.cos(lat) * np.sin(long)
    user_mcmf_pos_z = R_Moon * np.sin(lat)

    # 3x1 position vector
    user_mcmf_loc = np.array([user_mcmf_pos_x, user_mcmf_pos_y, user_mcmf_pos_z])
    return user_mcmf_loc

class EKF:
    def __init__(self, x0, P0) -> None:
        # define the process noise covariance
        self.Q = np.eye(3) * 0.1
        # define the measurement noise covariance
        self.R = np.eye(3) * 0.1
        # define the initial state estimate
        self.x_hat = x0
        # define the initial state estimate covariance
        self.P = P0

    def predict(self, rover, state, control):
        # predict the next state estimate
        A, B = rover.ss_model_dynamics(state, control)
        x_hat = rover.propagate_dynamics(state, control)
        P = A @ self.P @ A.T + self.Q
        return x_hat, P
    
    def update(self, rover, state, measurement):
        # update the state estimate
        #TODO: fix this to match the satellites
        # rover = Rover()
        H = np.eye(3)
        K = self.P @ H.T @ np.linalg.inv(H @ self.P @ H.T + self.R)
        x_hat = self.x_hat + K @ (measurement - state)
        P = (np.eye(3) - K @ H) @ self.P
        return x_hat, P
    



class MPC:
    def __init__(self, rover, horizon=10) -> None:
        self.rover = rover
        self.horizon = horizon

    def mpc_no_nav(
        x0: np.ndarray,
        A: np.ndarray,
        B: np.ndarray,
        P: np.ndarray,
        Q: np.ndarray,
        R: np.ndarray,
        W: np.ndarray,
        N: int,
        rx: float,
        ru: float,
    ) -> tuple[np.ndarray, np.ndarray, str]:
        """Solve the MPC problem starting at state `x0`."""
        n, m = Q.shape[0], R.shape[0]
        x_cvx = cvx.Variable((N + 1, n))
        u_cvx = cvx.Variable((N, m))

        # Construct and solve the MPC problem using CVXPY.
        # set up cost and constraints
        cost = 0
        constraints = [x_cvx[0] == x0]
        # iterate over time steps
        for k in range(N):
            cost = cost + cvx.quad_form(x_cvx[k], Q) + cvx.quad_form(u_cvx[k], R)
            # constraints on state and control values
            constraints += [cvx.norm(u_cvx[k]) <= ru]
            constraints += [cvx.norm(x_cvx[k]) <= rx]
            # dynamics constraints
            constraints += [x_cvx[k + 1] == A @ x_cvx[k] + B @ u_cvx[k]]

        # terminal constraint
        cost = cost + cvx.quad_form(x_cvx[N], P)
        constraints += [cvx.quad_form(x_cvx[N], W) <= 1]
        # solve the problem
        prob = cvx.Problem(cvx.Minimize(cost), constraints)
        prob.solve()
        x = x_cvx.value
        u = u_cvx.value
        status = prob.status
        return x, u, status
        # define the cost function

# def generate_deterministic_trajectory():
#     # A function that generates geometrically diverse trajectories
#     raise NotImplementedError


# def propagate_rover():
#     # A function that propagates the rover dynamics (discretized in time)
#     # Note, this might not be needed until after the midterm report
#     raise NotImplementedError