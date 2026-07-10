"""Reference-trajectory generators for the lunar-lander descent simulation (``ex11``).

Each generator returns a :class:`LanderTrajectory` sampled on the same time grid the C++
simulation uses (``N = round(T / dt) + 1`` epochs), in the **local East-North-Up frame** about
the DEM tile center, in meters (``U`` is height above the site datum). Feed the position array
into the lander app via ``set_reference_trajectory_enu`` (before ``sim.run()``) to make it the
truth trajectory:

    from pylupnt import lander_guidance as lg
    traj = lg.zem_zev_trajectory(r0, v0, rf, vf, T=300.0, dt=0.5)
    sim = pnt.Simulation("configs/lander_nav.yaml")
    app = sim.get_agent("Lander").get_application()
    app.set_reference_trajectory_enu(traj.r)   # (N, 3) ENU meters
    sim.run()

Three sources are provided:

* :func:`smoothstep_trajectory` -- a purely kinematic ``C1`` interpolation (the simulation's
  built-in default); no dynamics or thrust model.
* :func:`zem_zev_trajectory` -- **ZEM/ZEV feedback guidance**: the energy-optimal
  zero-effort-miss / zero-effort-velocity acceleration command, integrated forward.
* :func:`convex_descent_trajectory` -- **convex-optimization guidance**: a thrust-limited
  powered-descent solved as a second-order cone program with ``cvxpy`` (min-fuel) when it is
  installed, otherwise a ``scipy`` min-energy quadratic program on a coarse grid, interpolated
  to the fine grid.

The double-integrator model is ``r'' = a_thrust + g`` with ``g = (0, 0, -g_moon)``; the returned
``a`` is the total kinematic acceleration ``r''`` (so the commanded thrust is ``a - g``).
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

__all__ = [
    "LanderTrajectory",
    "MOON_SURFACE_GRAVITY",
    "smoothstep_trajectory",
    "zem_zev_trajectory",
    "convex_descent_trajectory",
]

# Mean lunar surface gravity [m/s^2] (GM_MOON / R_MOON^2).
MOON_SURFACE_GRAVITY = 1.625


class LanderTrajectory(NamedTuple):
    """A descent reference trajectory sampled on the simulation time grid.

    Attributes
    ----------
    t : (N,) ndarray
        Epoch times [s], ``t[k] = k * dt``.
    r : (N, 3) ndarray
        Position in local ENU about the tile center [m] (``U`` above the site datum).
    v : (N, 3) ndarray
        Velocity [m/s].
    a : (N, 3) ndarray
        Total kinematic acceleration ``r''`` [m/s^2]; commanded thrust is ``a - g``.
    """

    t: np.ndarray
    r: np.ndarray
    v: np.ndarray
    a: np.ndarray


def _grid(T: float, dt: float):
    n = max(2, int(round(T / dt)) + 1)
    return n, np.arange(n) * dt


def _as3(x) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size != 3:
        raise ValueError(f"expected a length-3 vector, got shape {np.shape(x)}")
    return x


def smoothstep_trajectory(r0, v0, rf, vf, T: float, dt: float) -> LanderTrajectory:
    """Kinematic smoothstep ``S(tau) = tau^2 (3 - 2 tau)`` interpolation from ``r0`` to ``rf``.

    Matches the simulation's built-in default. The path is ``C1`` with zero velocity and
    acceleration at both ends, so ``v0`` / ``vf`` are ignored (present only for a uniform
    signature). No dynamics or thrust limits are enforced.
    """
    r0, rf = _as3(r0), _as3(rf)
    n, t = _grid(T, dt)
    tau = np.clip(t / T, 0.0, 1.0)
    s = tau**2 * (3.0 - 2.0 * tau)
    sp = (6.0 * tau - 6.0 * tau**2) / T
    spp = (6.0 - 12.0 * tau) / T**2
    d = rf - r0
    r = r0[None, :] + np.outer(s, d)
    v = np.outer(sp, d)
    a = np.outer(spp, d)
    return LanderTrajectory(t, r, v, a)


def zem_zev_trajectory(
    r0,
    v0,
    rf,
    vf,
    T: float,
    dt: float,
    g: float = MOON_SURFACE_GRAVITY,
    tgo_min: float | None = None,
) -> LanderTrajectory:
    """Energy-optimal **ZEM/ZEV feedback guidance**, integrated forward on the time grid.

    At each step the time-to-go is ``tgo = T - t``; the zero-effort-miss and
    zero-effort-velocity (predicted terminal errors under gravity-only coasting) are

    .. code::

        ZEM = rf - (r + v*tgo + 0.5*g*tgo^2)
        ZEV = vf - (v + g*tgo)

    and the commanded thrust acceleration is ``a_thrust = (6/tgo^2) ZEM - (2/tgo) ZEV`` (the
    minimum-energy solution of the double integrator). The vehicle acceleration ``r'' =
    a_thrust + g`` is integrated with a semi-implicit step. ``tgo`` is floored at ``tgo_min``
    (default ``max(2*dt, 2% of T)``) so the terminal gains stay finite.

    Parameters
    ----------
    r0, v0, rf, vf : length-3
        Initial / final ENU position [m] and velocity [m/s].
    g : float
        Lunar gravity magnitude [m/s^2] (downward, ``-U``).
    """
    r0, v0, rf, vf = _as3(r0), _as3(v0), _as3(rf), _as3(vf)
    n, t = _grid(T, dt)
    g_vec = np.array([0.0, 0.0, -g])
    if tgo_min is None:
        tgo_min = max(2.0 * dt, 0.02 * T)

    r = np.zeros((n, 3))
    v = np.zeros((n, 3))
    a = np.zeros((n, 3))
    r[0], v[0] = r0, v0
    for k in range(n):
        tgo = max(T - t[k], tgo_min)
        zem = rf - (r[k] + v[k] * tgo + 0.5 * g_vec * tgo**2)
        zev = vf - (v[k] + g_vec * tgo)
        a_thrust = (6.0 / tgo**2) * zem - (2.0 / tgo) * zev
        a_tot = a_thrust + g_vec
        a[k] = a_tot
        if k + 1 < n:
            v[k + 1] = v[k] + a_tot * dt
            r[k + 1] = r[k] + v[k] * dt + 0.5 * a_tot * dt**2
    return LanderTrajectory(t, r, v, a)


def _propagate(a_thrust: np.ndarray, r0, v0, g: float, h: float):
    """Forward-propagate the double integrator ``r'' = a_thrust + g`` (ZOH thrust, step ``h``)."""
    m = a_thrust.shape[0]
    g_vec = np.array([0.0, 0.0, -g])
    r = np.zeros((m, 3))
    v = np.zeros((m, 3))
    r[0], v[0] = r0, v0
    for k in range(m - 1):
        acc = a_thrust[k] + g_vec
        v[k + 1] = v[k] + acc * h
        r[k + 1] = r[k] + v[k] * h + 0.5 * acc * h**2
    return r, v


def _convex_cvxpy(r0, v0, rf, vf, T, g, a_max, n_nodes):
    import cvxpy as cp

    m = n_nodes
    h = T / (m - 1)
    g_vec = np.array([0.0, 0.0, -g])
    A = cp.Variable((m, 3))  # thrust acceleration at each node
    r = [np.asarray(r0)]
    v = [np.asarray(v0)]
    for k in range(m - 1):
        acc = A[k, :] + g_vec
        v.append(v[k] + acc * h)
        r.append(r[k] + v[k] * h + 0.5 * acc * h**2)
    cons = [r[-1] == rf, v[-1] == vf]
    cons += [cp.norm(A[k, :]) <= a_max for k in range(m)]  # SOC thrust-magnitude limit
    obj = cp.Minimize(cp.sum([cp.norm(A[k, :]) for k in range(m)]) * h)  # min-fuel proxy
    prob = cp.Problem(obj, cons)
    prob.solve()
    if A.value is None:
        raise RuntimeError(f"cvxpy solve failed (status={prob.status})")
    return np.asarray(A.value)


def _convex_scipy(r0, v0, rf, vf, T, g, a_max, n_nodes):
    from scipy.optimize import minimize

    m = n_nodes
    h = T / (m - 1)

    def terminal_err(x):
        a = x.reshape(m, 3)
        r, v = _propagate(a, r0, v0, g, h)
        return np.concatenate([r[-1] - rf, v[-1] - vf])

    def energy(x):
        return 0.5 * float(np.sum(x**2))

    def energy_grad(x):
        return x

    # Thrust-magnitude cone as smooth inequalities a_max^2 - |a_k|^2 >= 0.
    def thrust_ineq(x):
        a = x.reshape(m, 3)
        return a_max**2 - np.sum(a**2, axis=1)

    # Warm-start from the (unconstrained) ZEM/ZEV thrust profile sampled at the coarse nodes.
    _, tc = _grid(T, h)
    zz = zem_zev_trajectory(r0, v0, rf, vf, T, h, g=g)
    x0 = (zz.a - np.array([0.0, 0.0, -g]))[:m].reshape(-1)

    cons = [
        {"type": "eq", "fun": terminal_err},
        {"type": "ineq", "fun": thrust_ineq},
    ]
    res = minimize(
        energy,
        x0,
        jac=energy_grad,
        constraints=cons,
        method="SLSQP",
        options={"maxiter": 400, "ftol": 1e-6},
    )
    if not res.success and np.linalg.norm(terminal_err(res.x)) > 1.0:
        raise RuntimeError(f"scipy SLSQP solve failed: {res.message}")
    return res.x.reshape(m, 3)


def convex_descent_trajectory(
    r0,
    v0,
    rf,
    vf,
    T: float,
    dt: float,
    g: float = MOON_SURFACE_GRAVITY,
    a_max: float = 5.0,
    n_nodes: int = 40,
    solver: str = "auto",
) -> LanderTrajectory:
    """Thrust-limited **convex-optimization** powered descent, on the simulation time grid.

    Solves a double-integrator descent with a thrust-magnitude limit ``|a_thrust| <= a_max`` on
    a coarse grid of ``n_nodes`` nodes, then cubic-interpolates the position to the fine grid.
    With ``cvxpy`` installed the problem is the standard **min-fuel second-order cone program**
    (minimize ``sum |a_thrust| dt`` s.t. dynamics, boundary conditions, and the thrust cone);
    otherwise a ``scipy`` SLSQP **min-energy** solve (minimize ``sum |a_thrust|^2``) with the
    same constraints is used, warm-started from the ZEM/ZEV profile.

    Parameters
    ----------
    a_max : float
        Maximum thrust acceleration magnitude [m/s^2].
    n_nodes : int
        Number of coarse optimization nodes (kept small so the solve is fast; the trajectory is
        smooth, so interpolation to the fine grid is accurate).
    solver : {"auto", "cvxpy", "scipy"}
        Which backend to use. ``"auto"`` tries ``cvxpy`` and falls back to ``scipy``.
    """
    from scipy.interpolate import CubicSpline

    r0, v0, rf, vf = _as3(r0), _as3(v0), _as3(rf), _as3(vf)
    m = max(4, int(n_nodes))
    h = T / (m - 1)

    if solver == "cvxpy":
        a_coarse = _convex_cvxpy(r0, v0, rf, vf, T, g, a_max, m)
    elif solver == "scipy":
        a_coarse = _convex_scipy(r0, v0, rf, vf, T, g, a_max, m)
    else:  # auto
        try:
            a_coarse = _convex_cvxpy(r0, v0, rf, vf, T, g, a_max, m)
        except Exception:
            a_coarse = _convex_scipy(r0, v0, rf, vf, T, g, a_max, m)

    r_coarse, _ = _propagate(a_coarse, r0, v0, g, h)
    t_coarse = np.arange(m) * h

    n, t = _grid(T, dt)
    spl = [CubicSpline(t_coarse, r_coarse[:, i]) for i in range(3)]
    r = np.column_stack([s(t) for s in spl])
    v = np.column_stack([s(t, 1) for s in spl])
    a = np.column_stack([s(t, 2) for s in spl])
    return LanderTrajectory(t, r, v, a)
