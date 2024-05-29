import os
import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import pylupnt as pnt
import scipy as sc


class Rover_Agent:
    def __init__(self):
        # constants
        self.g = 1.62       # m/s^2
        self.m = 487        # kg (Endurance report states 487 fully margined)
        self.dt = 0.01      # s
        self.N = 50         # grid size

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
    


# class EKF:
#     def __init__(self) -> None:
#         # define the process noise covariance
#         self.Q = jnp.eye(3) * 0.1
#         # define the measurement noise covariance
#         self.R = jnp.eye(3) * 0.1
#         # define the initial state estimate
#         self.x_hat = jnp.array([0, 0, 0])
#         # define the initial state estimate covariance
#         self.P = jnp.eye(3) * 0.1

#     def predict(self, state, control):
#         # predict the next state estimate
#         rover = Rover()
#         A, B = rover.ss_model_dynamics(state, control)
#         x_hat = rover.propagate_dynamics(state, control)
#         P = A @ self.P @ A.T + self.Q
#         return x_hat, P
    
#     def update(self, state, measurement):
#         # update the state estimate
#         rover = Rover()
#         H = jnp.eye(3)
#         K = self.P @ H.T @ jnp.linalg.inv(H @ self.P @ H.T + self.R)
#         x_hat = self.x_hat + K @ (measurement - state)
#         P = (jnp.eye(3) - K @ H) @ self.P
#         return x_hat, P
    






# def generate_deterministic_trajectory():
#     # A function that generates geometrically diverse trajectories
#     raise NotImplementedError


# def propagate_rover():
#     # A function that propagates the rover dynamics (discretized in time)
#     # Note, this might not be needed until after the midterm report
#     raise NotImplementedError