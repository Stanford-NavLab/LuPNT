"""Unit tests for the pure-Python lunar-lander reference-trajectory generators in
``pylupnt.lander_guidance`` (ex11): the ``_as3``/``_grid`` helpers, the kinematic
smoothstep path, ZEM/ZEV feedback guidance, forward propagation, and the convex
(cvxpy / scipy) powered-descent guidance.
"""

import sys
import types

import numpy as np
import pytest

from pylupnt import lander_guidance as lg


# --------------------------------------------------------------------------- constants / helpers
def test_moon_surface_gravity_value():
    assert lg.MOON_SURFACE_GRAVITY == pytest.approx(1.625)


def test_as3_accepts_various_shapes():
    np.testing.assert_array_equal(lg._as3([1, 2, 3]), [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(lg._as3(np.array([[4.0], [5.0], [6.0]])), [4.0, 5.0, 6.0])


def test_as3_rejects_wrong_length():
    with pytest.raises(ValueError, match="length-3"):
        lg._as3([1.0, 2.0])


def test_grid_minimum_two_points_and_spacing():
    n, t = lg._grid(10.0, 2.0)
    assert n == 6  # round(10/2)+1
    np.testing.assert_allclose(t, [0, 2, 4, 6, 8, 10])
    # Degenerate T -> still at least 2 nodes.
    n2, t2 = lg._grid(0.0, 1.0)
    assert n2 == 2


# --------------------------------------------------------------------------- smoothstep_trajectory
def test_smoothstep_endpoints_and_zero_end_velocity():
    r0 = np.array([0.0, 0.0, 100.0])
    rf = np.array([50.0, -20.0, 0.0])
    traj = lg.smoothstep_trajectory(r0, np.zeros(3), rf, np.zeros(3), T=60.0, dt=1.0)

    assert traj.r.shape == traj.v.shape == traj.a.shape == (61, 3)
    assert traj.t.shape == (61,)
    np.testing.assert_allclose(traj.r[0], r0, atol=1e-9)
    np.testing.assert_allclose(traj.r[-1], rf, atol=1e-9)
    # Smoothstep has zero velocity at both ends.
    np.testing.assert_allclose(traj.v[0], np.zeros(3), atol=1e-9)
    np.testing.assert_allclose(traj.v[-1], np.zeros(3), atol=1e-9)


def test_smoothstep_is_monotone_along_each_axis():
    r0 = np.zeros(3)
    rf = np.array([10.0, 0.0, 0.0])
    traj = lg.smoothstep_trajectory(r0, None, rf, None, T=30.0, dt=0.5)
    x = traj.r[:, 0]
    assert np.all(np.diff(x) >= -1e-9)  # non-decreasing
    assert x[0] == pytest.approx(0.0)
    assert x[-1] == pytest.approx(10.0)


# --------------------------------------------------------------------------- zem_zev_trajectory
def test_zem_zev_hits_boundary_conditions():
    r0 = np.array([0.0, 0.0, 200.0])
    v0 = np.array([5.0, 0.0, -10.0])
    rf = np.array([0.0, 0.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])
    traj = lg.zem_zev_trajectory(r0, v0, rf, vf, T=120.0, dt=0.2)

    np.testing.assert_allclose(traj.r[0], r0, atol=1e-9)
    np.testing.assert_allclose(traj.v[0], v0, atol=1e-9)
    # Terminal state should be close to the target (integrated feedback law).
    np.testing.assert_allclose(traj.r[-1], rf, atol=1.0)
    np.testing.assert_allclose(traj.v[-1], vf, atol=0.5)
    assert traj.a.shape == (601, 3)


def test_zem_zev_explicit_tgo_min_and_custom_gravity():
    r0 = np.array([10.0, 5.0, 80.0])
    v0 = np.zeros(3)
    rf = np.zeros(3)
    vf = np.zeros(3)
    # Explicit tgo_min branch + non-default gravity magnitude.
    traj = lg.zem_zev_trajectory(r0, v0, rf, vf, T=40.0, dt=0.5, g=1.0, tgo_min=5.0)
    assert np.all(np.isfinite(traj.a))
    np.testing.assert_allclose(traj.r[-1], rf, atol=1.0)


# --------------------------------------------------------------------------- _propagate
def test_propagate_matches_analytic_gravity_freefall():
    # Zero thrust: pure free-fall in -z. r(t) = r0 - 0.5 g t^2.
    m = 5
    h = 1.0
    g = 1.625
    a_thrust = np.zeros((m, 3))
    r0 = np.array([0.0, 0.0, 100.0])
    v0 = np.zeros(3)
    r, v = lg._propagate(a_thrust, r0, v0, g, h)
    assert r.shape == (m, 3) and v.shape == (m, 3)
    # Semi-implicit (symplectic Euler) integrator with ZOH acceleration.
    t = np.arange(m) * h
    # Velocity is exact for constant acceleration: v_z = -g t.
    np.testing.assert_allclose(v[:, 2], -g * t, atol=1e-9)
    # Falling (monotonically decreasing height).
    assert np.all(np.diff(r[:, 2]) <= 0)


# --------------------------------------------------------------------------- convex_descent (scipy)
def test_convex_descent_scipy_hits_boundary_conditions():
    r0 = np.array([0.0, 0.0, 100.0])
    v0 = np.array([0.0, 0.0, -5.0])
    rf = np.zeros(3)
    vf = np.zeros(3)
    traj = lg.convex_descent_trajectory(
        r0, v0, rf, vf, T=60.0, dt=1.0, g=1.625, a_max=5.0, n_nodes=6, solver="scipy"
    )
    assert traj.r.shape == (61, 3)
    np.testing.assert_allclose(traj.r[0], r0, atol=2.0)
    np.testing.assert_allclose(traj.r[-1], rf, atol=2.0)
    assert np.all(np.isfinite(traj.a))


def test_convex_descent_auto_falls_back_to_scipy_without_cvxpy():
    # cvxpy is not installed -> the "auto" branch catches the ImportError and
    # falls back to the scipy solve.
    if "cvxpy" in sys.modules:
        pytest.skip("cvxpy present; auto path would use it")
    r0 = np.array([0.0, 0.0, 80.0])
    v0 = np.zeros(3)
    rf = np.zeros(3)
    vf = np.zeros(3)
    traj = lg.convex_descent_trajectory(r0, v0, rf, vf, T=50.0, dt=1.0, n_nodes=6, solver="auto")
    assert traj.r.shape[0] == 51
    np.testing.assert_allclose(traj.r[-1], rf, atol=2.0)


# --------------------------------------------------------------------------- convex_descent (cvxpy)
class _CpExpr:
    """Minimal stand-in for a cvxpy expression that swallows all arithmetic."""

    __array_ufunc__ = None  # let numpy defer to our reflected operators

    def __add__(self, other):
        return self

    __radd__ = __add__
    __sub__ = __add__
    __rsub__ = __add__

    def __mul__(self, other):
        return self

    __rmul__ = __mul__

    def __eq__(self, other):
        return ("eq", other)

    def __le__(self, other):
        return ("le", other)

    def __hash__(self):
        return id(self)


class _CpVar(_CpExpr):
    def __init__(self, shape):
        self.shape = shape
        self.value = np.zeros(shape)  # non-None so the solve check passes

    def __getitem__(self, idx):
        return _CpExpr()


def _make_fake_cvxpy(leave_value=True):
    mod = types.ModuleType("cvxpy")
    holder = {}

    def Variable(shape):
        v = _CpVar(shape)
        holder["var"] = v
        return v

    class Problem:
        def __init__(self, obj, cons):
            self.status = "optimal"

        def solve(self):
            if not leave_value:
                holder["var"].value = None
            return 0.0

    mod.Variable = Variable
    mod.norm = lambda x: _CpExpr()
    mod.sum = lambda seq: _CpExpr()
    mod.Minimize = lambda x: x
    mod.Problem = Problem
    return mod


def test_convex_descent_cvxpy_path_with_fake_solver(monkeypatch):
    monkeypatch.setitem(sys.modules, "cvxpy", _make_fake_cvxpy(leave_value=True))
    r0 = np.array([0.0, 0.0, 90.0])
    v0 = np.zeros(3)
    rf = np.zeros(3)
    vf = np.zeros(3)
    traj = lg.convex_descent_trajectory(r0, v0, rf, vf, T=50.0, dt=1.0, n_nodes=6, solver="cvxpy")
    assert traj.r.shape == (51, 3)
    assert np.all(np.isfinite(traj.a))


def test_convex_cvxpy_raises_when_solver_leaves_no_value(monkeypatch):
    monkeypatch.setitem(sys.modules, "cvxpy", _make_fake_cvxpy(leave_value=False))
    with pytest.raises(RuntimeError, match="cvxpy solve failed"):
        lg._convex_cvxpy(
            np.zeros(3),
            np.zeros(3),
            np.zeros(3),
            np.zeros(3),
            T=50.0,
            g=1.625,
            a_max=5.0,
            n_nodes=6,
        )
