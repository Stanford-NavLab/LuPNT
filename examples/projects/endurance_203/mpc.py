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

    # Initialize trajectory
    # if s_init is None or u_init is None:
    #     s = np.zeros((N + 1, n))
    #     u = np.zeros((N, m))
    #     # u = np.random.randn(N,m)
    #     s[0] = s0
    #     # print('In here (initialization)')
    #     for k in range(N):
    #         s[k + 1] = f(s[k], u[k])
    # else:
    #     s = np.copy(s_init)
    #     u = np.copy(u_init)
    
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

def run_mpc(s0, s_goal, N, P, Q, R, T, N_scp, dt, v_bound, omega_bound):
    n = 3 # state dimension
    m = 2 # control dimension
    eps = 1e-3 # SCP convergence tolerance

    # N = 5 # MPC horizon length
    # # P = 1e2 * np.eye(n)  # terminal state cost matrix
    # # Q = np.eye(n)
    # P = np.diag([10, 10, 1])*1e2
    # Q = np.diag([1, 1, 1])
    # R = 1e-2 * np.eye(m)  # control cost matrix
    # s0 = np.array([2.0, 3.0, -np.pi/2]) # initial state (should be previous waypoint)
    # s_goal = np.array([7.0, 5.0, -np.pi/4]) # desired final state (should be next waypoint)
    # T = 30 # total simulation time
    # N_scp = 10 # maximum number of SCP iterations
    # dt = 0.1 # time step

    f = dynamics
    s_mpc = np.zeros((T, N + 1, n))
    u_mpc = np.zeros((T, N, m))
    s = np.copy(s0)
    total_time = time()
    total_control_cost = 0.0
    s_init = None
    u_init = None

    for t in tqdm(range(T)):
        s_mpc[t], u_mpc[t] = solve_mpc(f, s, s_goal, N, P, Q, R, eps, N_scp, dt, v_bound, omega_bound, s_init, u_init)
        s = f(s, u_mpc[t, 0])
        # print(f'Applied u = {u_mpc[t, 0]}')
        
        total_control_cost += u_mpc[t, 0].T @ R @ u_mpc[t, 0]

        # Use this solution to warm-start the next iteration
        u_init = np.concatenate([u_mpc[t, 1:], u_mpc[t, -1:]])
        s_init = np.concatenate(
            [s_mpc[t, 1:], f(s_mpc[t, -1], u_mpc[t, -1]).reshape([1, -1])]
        )
    total_time = time() - total_time
    print("Total elapsed time:", total_time, "seconds")
    print("Total control cost:", total_control_cost)
    # print(s_mpc[:, 0, 1])
    return s_mpc, u_mpc

def plot_mpc(s0, s_goal, s_mpc, u_mpc, N, T, N_scp, n_waypt, dt):
    fig, ax = plt.subplots(2, 2, dpi=150, figsize=(15, 10))
    # fig.suptitle("$N = {}$, ".format(N) + r"$N_\mathrm{SCP} = " + "{}$".format(N_scp))
    fig.suptitle('State and Control Input Over Time')

    for t in range(T):
        ax[0, 0].plot(s_mpc[t, :, 1], s_mpc[t, :, 0], "--*", color="k")
    ax[0, 0].plot(s_mpc[:-1, 0, 1], s_mpc[:-1, 0, 0], "-o")
    ax[0, 0].set_xlabel(r"$y(t)$")
    ax[0, 0].set_ylabel(r"$x(t)$")
    ax[0, 0].axis("equal")
    ax[0, 0].scatter(s0[1], s0[0],color='r',zorder=3)
    ax[0, 0].scatter(s_goal[1], s_goal[0],color='r',zorder=3)

    t_vec = np.linspace(0, T*(n_waypt - 1)*dt, len(s_mpc[:, 0, 2]))
    ax[0, 1].plot(t_vec, s_mpc[:, 0, 2]*180/np.pi, "-o") #, label=r"$\theta(t)$")
    ax[0, 1].set_xlabel(r"$t$ [hr]")
    ax[0, 1].set_ylabel(r"$\theta(t)$ [deg]")
    # ax[0, 1].legend()
    # ax[0, 1].axhline(s_goal[2]*180/np.pi)

    ax[1, 0].plot(t_vec, u_mpc[:, 0, 0], "-o") #, label=r"$u_1(t)$")
    ax[1, 0].set_xlabel(r"$t$ [hr]")
    ax[1, 0].set_ylabel(r"$v(t)$ [km/hr]")
    # ax[1, 0].legend()
    ax[1, 0].set_ylim([-0.1, 5.1])

    ax[1, 1].plot(t_vec, u_mpc[:, 0, 1]*180/np.pi, "-o") #, label=r"$u_2(t)$")
    ax[1, 1].set_xlabel(r"$t$ [hr]")
    ax[1, 1].set_ylabel(r"$\omega(t)$ [deg/hr]")
    # ax[1, 1].legend()
    ax[1, 1].set_ylim(np.array([-3.1, 3.1])*180/np.pi)

    # suffix = "_N={}_Nscp={}".format(N, N_scp)
    plt.tight_layout()
    # plt.savefig("soln_obstacle_avoidance" + suffix + ".png", bbox_inches="tight")
    plt.savefig("mpc_state_control_over_time.png", bbox_inches="tight")

    # plt.show()
