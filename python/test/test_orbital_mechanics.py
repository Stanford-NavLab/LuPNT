"""Unit tests for `pylupnt` orbital-mechanics bindings:
`get_orbital_period`, `propagate_coe` (Keplerian propagation),
`tle_to_classical`, and the classical <-> Cartesian round trip.
"""

import os

import numpy as np
import pylupnt as pnt
import pytest

# A circular-ish lunar orbit: [a, e, i, raan, argp, M] (SI + radians).
COE0 = np.array([6541.4e3, 0.05, np.radians(56.2), 0.0, np.radians(90.0), 0.0])
GM = pnt.GM_MOON


# ---------------------------------------------------------------------------
# get_orbital_period
# ---------------------------------------------------------------------------


def test_orbital_period_matches_kepler_third_law():
    a = COE0[0]
    period = float(pnt.get_orbital_period(a, GM))
    expected = 2 * np.pi * np.sqrt(a**3 / GM)
    assert period == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# propagate_coe
# ---------------------------------------------------------------------------


def test_propagate_coe_conserves_shape_elements():
    dt = 1234.0
    coe = np.asarray(pnt.propagate_coe(COE0, GM, dt)).reshape(-1)
    # Two-body propagation only advances the mean anomaly; a, e, i, RAAN, argp fixed.
    np.testing.assert_allclose(coe[:5], COE0[:5], atol=1e-6)
    assert coe[5] != COE0[5]


def test_propagate_coe_full_period_returns_to_start():
    period = float(pnt.get_orbital_period(COE0[0], GM))
    coe = np.asarray(pnt.propagate_coe(COE0, GM, period)).reshape(-1)
    # After one period the mean anomaly returns to its start (mod 2*pi).
    dm = (coe[5] - COE0[5] + np.pi) % (2 * np.pi) - np.pi
    assert dm == pytest.approx(0.0, abs=1e-6)


def test_propagate_coe_mean_anomaly_rate():
    dt = 600.0
    n = np.sqrt(GM / COE0[0] ** 3)  # mean motion [rad/s]
    coe = np.asarray(pnt.propagate_coe(COE0, GM, dt)).reshape(-1)
    dm = (coe[5] - COE0[5]) % (2 * np.pi)
    assert dm == pytest.approx((n * dt) % (2 * np.pi), abs=1e-6)


# ---------------------------------------------------------------------------
# classical <-> Cartesian
# ---------------------------------------------------------------------------


def test_classical_cartesian_round_trip():
    rv = np.asarray(pnt.classical_to_cart(COE0, GM))
    assert rv.shape == (6,)
    coe = np.asarray(pnt.cart_to_classical(rv, GM)).reshape(-1)
    np.testing.assert_allclose(coe[:5], COE0[:5], atol=1e-3)


# ---------------------------------------------------------------------------
# TLE -> classical elements
# ---------------------------------------------------------------------------


def _gps_tle_file():
    return os.path.join(pnt.get_data_path(), "tle", "2025_01_01", "gps_2025_01_01.txt")


def test_tle_to_classical_gps():
    tle_file = _gps_tle_file()
    if not os.path.isfile(tle_file):
        pytest.skip(f"TLE test data not found at {tle_file}")

    lines = [ln.strip() for ln in open(tle_file) if ln.strip()]
    tle = pnt.TLE.from_lines(lines[0], lines[1], lines[2])
    t_tai = tle.epoch_tai  # property [TAI seconds]

    coe = np.asarray(pnt.tle_to_classical(tle, t_tai)).reshape(-1)
    a_km = coe[0] / 1e3
    inc_deg = np.degrees(coe[2])
    assert 26_000 < a_km < 27_500  # GPS semi-major axis ~26,560 km
    assert 0.0 <= coe[1] < 0.05  # near-circular
    assert 50 < inc_deg < 60  # GPS inclination ~55 deg
