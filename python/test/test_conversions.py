"""Unit tests for the `pylupnt` coordinate/state conversion bindings
(`python/bindings/py_conversions.cc`): anomaly conversions, attitude/
quaternion helpers, topocentric (ENU/AER) and geodetic (LLA/stereographic)
transforms, and the classical/quasi-nonsingular/equinoctial/Delaunay orbital
element sets.

These exercise the same `lupnt` free functions the C++ library uses, through
their Python bindings. Because every listed transform has a documented inverse,
the tests are built around round-trip invariants (`f_inv(f(x)) == x`) plus a
handful of closed-form spot checks (Kepler's equation, unit vectors, det = +1).
"""

import numpy as np
import pylupnt as pnt
import pytest

ABS_TOL = 1e-6


def _f(x):
    """Coerce an autodiff-backed Vec/Mat return into a plain float ndarray."""
    return np.asarray(x, dtype=float)


# ---------------------------------------------------------------------------
# Anomaly conversions
# ---------------------------------------------------------------------------


def test_orbital_period():
    # Circular orbit period T = 2*pi*sqrt(a^3/GM)
    a, GM = 5740e3, pnt.GM_MOON
    expected = pnt.TWO_PI * np.sqrt(a**3 / GM)
    assert float(pnt.get_orbital_period(a, GM)) == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("e", [0.0, 0.1, 0.58, 0.9])
@pytest.mark.parametrize("M", [0.0, 0.3, 1.5, 3.0])
def test_anomaly_round_trips(e, M):
    # M -> E solves Kepler's equation E - e*sin(E) = M
    E = float(pnt.mean_to_eccentric_anomaly(M, e))
    assert E - e * np.sin(E) == pytest.approx(M, abs=ABS_TOL)

    # E -> M is the forward (closed-form) direction and inverts M -> E
    assert float(pnt.eccentric_to_mean_anomaly(E, e)) == pytest.approx(M, abs=ABS_TOL)

    # E <-> nu (true anomaly) is a bijection on the branch near M
    nu = float(pnt.eccentric_to_true_anomaly(E, e))
    assert float(pnt.true_to_eccentric_anomaly(nu, e)) == pytest.approx(E, abs=ABS_TOL)

    # The M <-> nu shortcuts compose the two chained conversions above
    assert float(pnt.mean_to_true_anomaly(M, e)) == pytest.approx(nu, abs=ABS_TOL)
    assert float(pnt.true_to_mean_anomaly(nu, e)) == pytest.approx(M, abs=ABS_TOL)


def test_anomaly_circular_identity():
    # For a circular orbit (e = 0) mean = eccentric = true anomaly
    for M in (0.0, 0.7, 2.0):
        assert float(pnt.mean_to_eccentric_anomaly(M, 0.0)) == pytest.approx(M, abs=ABS_TOL)
        assert float(pnt.mean_to_true_anomaly(M, 0.0)) == pytest.approx(M, abs=ABS_TOL)


# ---------------------------------------------------------------------------
# Attitude / quaternion conversions
# ---------------------------------------------------------------------------


def test_quaternion_conventions():
    q = np.array([1.0, 2.0, 3.0, 4.0])

    # normalize_quat yields a unit quaternion
    qn = _f(pnt.normalize_quat(q))
    assert np.linalg.norm(qn) == pytest.approx(1.0, abs=ABS_TOL)

    # scalar_first_to_last / scalar_last_to_first permute the components between
    # the two ordering conventions (both normalize their input) and, being
    # opposite reorderings, invert one another
    q_last = _f(pnt.scalar_first_to_last(qn))
    assert np.linalg.norm(q_last) == pytest.approx(1.0, abs=ABS_TOL)
    assert sorted(np.abs(q_last)) == pytest.approx(sorted(np.abs(qn)), abs=ABS_TOL)
    np.testing.assert_allclose(_f(pnt.scalar_last_to_first(q_last)), qn, atol=ABS_TOL)
    np.testing.assert_allclose(
        _f(pnt.scalar_first_to_last(_f(pnt.scalar_last_to_first(qn)))), qn, atol=ABS_TOL
    )


def test_quat_rot_round_trip():
    qn = _f(pnt.normalize_quat(np.array([1.0, 2.0, 3.0, 4.0])))

    R = pnt.quat_to_rot(qn)
    # A rotation matrix is orthonormal with determinant +1
    np.testing.assert_allclose(R.T @ R, np.eye(3), atol=ABS_TOL)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=ABS_TOL)

    # rot_to_quat recovers the quaternion up to an overall sign (q and -q are
    # the same rotation)
    q2 = _f(pnt.rot_to_quat(R))
    np.testing.assert_allclose(np.abs(q2), np.abs(qn), atol=ABS_TOL)


def test_roll_pitch_yaw_round_trip():
    rpy = np.array([0.1, -0.2, 0.3])
    R = pnt.roll_pitch_yaw_to_rot(rpy)
    np.testing.assert_allclose(R.T @ R, np.eye(3), atol=ABS_TOL)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=ABS_TOL)
    np.testing.assert_allclose(_f(pnt.rot_to_roll_pitch_yaw(R)), rpy, atol=ABS_TOL)


# ---------------------------------------------------------------------------
# Topocentric conversions (ENU <-> AER <-> Cartesian)
# ---------------------------------------------------------------------------


def test_enu_aer_round_trips():
    xyz_ref = np.array([pnt.R_MOON, 0.0, 0.0])
    xyz = xyz_ref + np.array([100.0, 200.0, 300.0])

    # Cartesian -> ENU -> Cartesian
    enu = pnt.cart_to_east_north_up(xyz, xyz_ref)
    np.testing.assert_allclose(_f(pnt.east_north_up_to_cart(enu, xyz_ref)), xyz, atol=ABS_TOL)

    # ENU -> AER -> ENU; range is the ENU vector magnitude
    aer = pnt.east_north_up_to_az_el_range(enu)
    assert _f(aer)[2] == pytest.approx(np.linalg.norm(_f(enu)), abs=ABS_TOL)
    np.testing.assert_allclose(_f(pnt.az_el_range_to_east_north_up(aer)), _f(enu), atol=ABS_TOL)

    # The compound AER <-> Cartesian helpers agree with going through ENU
    np.testing.assert_allclose(_f(pnt.cart_to_az_el_range(xyz, xyz_ref)), _f(aer), atol=ABS_TOL)
    np.testing.assert_allclose(_f(pnt.az_el_range_to_cart(aer, xyz_ref)), xyz, atol=ABS_TOL)


# ---------------------------------------------------------------------------
# Geodetic conversions (LLA / stereographic <-> Cartesian)
# ---------------------------------------------------------------------------


def test_lat_lon_alt_round_trips():
    R = pnt.R_MOON
    lla = np.array([0.3, 0.5, 1000.0])  # [rad, rad, m]

    xyz = pnt.lat_lon_alt_to_cart(lla, R)
    # A point at altitude h sits at radius R + h from the center
    assert np.linalg.norm(_f(xyz)) == pytest.approx(R + lla[2], abs=1e-3)
    np.testing.assert_allclose(_f(pnt.cart_to_lat_lon_alt(xyz, R)), lla, atol=ABS_TOL)


def test_stereographic_round_trips():
    R = pnt.R_MOON
    lla = np.array([0.3, 0.5, 1000.0])

    xya = pnt.lat_lon_alt_to_stereographic(lla, R)
    np.testing.assert_allclose(_f(pnt.stereographic_to_lat_lon_alt(xya, R)), lla, atol=ABS_TOL)

    # stereographic <-> Cartesian is consistent with the LLA <-> Cartesian path
    xyz = pnt.stereographic_to_cart(xya, R)
    np.testing.assert_allclose(_f(pnt.cart_to_stereographic(xyz, R)), _f(xya), atol=ABS_TOL)
    np.testing.assert_allclose(_f(pnt.lat_lon_alt_to_cart(lla, R)), _f(xyz), atol=1e-3)


# ---------------------------------------------------------------------------
# Orbital element-set conversions
# ---------------------------------------------------------------------------


# [a (m), e, i (rad), RAAN (rad), argp (rad), M (rad)]
_COE = np.array([5740e3, 0.58, 0.9, 0.1, 1.5, 0.3])


def test_classical_cartesian_round_trip():
    GM = pnt.GM_MOON
    rv = pnt.classical_to_cart(_COE, GM)
    coe_back = _f(pnt.cart_to_classical(rv, GM))
    np.testing.assert_allclose(coe_back, _COE, atol=ABS_TOL)

    # Sanity: the Cartesian radius stays within [perigee, apogee]
    a, e = _COE[0], _COE[1]
    r = np.linalg.norm(_f(rv)[:3])
    assert a * (1 - e) - 1.0 <= r <= a * (1 + e) + 1.0


@pytest.mark.parametrize(
    "to_fn,from_fn",
    [
        ("classical_to_quasi_nonsingular", "quasi_nonsingular_to_classical"),
        ("classical_to_equinoctial", "equinoctial_to_classical"),
        ("classical_to_delaunay", "delaunay_to_classical"),
    ],
)
def test_element_set_round_trips(to_fn, from_fn):
    GM = pnt.GM_MOON
    other = getattr(pnt, to_fn)(_COE, GM)
    coe_back = _f(getattr(pnt, from_fn)(other, GM))
    np.testing.assert_allclose(coe_back, _COE, atol=ABS_TOL)


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
