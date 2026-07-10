"""Unit tests for `pylupnt` math / unit helper bindings:
`skew`, `rot_x`/`rot_y`/`rot_z`, `wrap_to_pi`/`wrap_to_two_pi`,
`decibel_to_decimal`/`decimal_to_decibel`,
`degrees_to_deg_min_sec`/`deg_min_sec_to_degrees`.
"""

import numpy as np
import pylupnt as pnt
import pytest


# ---------------------------------------------------------------------------
# skew
# ---------------------------------------------------------------------------


def test_skew_reproduces_cross_product():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([-4.0, 0.5, 2.0])
    S = np.asarray(pnt.skew(a))
    assert S.shape == (3, 3)
    # skew(a) is antisymmetric and skew(a) @ b == a x b
    np.testing.assert_allclose(S, -S.T, atol=1e-12)
    np.testing.assert_allclose(S @ b, np.cross(a, b), atol=1e-12)


# ---------------------------------------------------------------------------
# Elementary rotation matrices
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("rot", ["rot_x", "rot_y", "rot_z"])
def test_rotation_matrices_are_orthonormal(rot):
    R = np.asarray(getattr(pnt, rot)(0.7))
    assert R.shape == (3, 3)
    np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)


def test_rot_z_known_angle():
    R = np.asarray(pnt.rot_z(np.pi / 2))
    # rot_z(90 deg) maps x-hat to +/- y-hat (sign is convention dependent).
    mapped = R @ np.array([1.0, 0.0, 0.0])
    assert abs(mapped[2]) < 1e-12
    assert abs(abs(mapped[1]) - 1.0) < 1e-12


def test_rotation_matrices_compose_back_to_identity():
    R = np.asarray(pnt.rot_x(0.3))
    Rinv = np.asarray(pnt.rot_x(-0.3))
    np.testing.assert_allclose(R @ Rinv, np.eye(3), atol=1e-12)


# ---------------------------------------------------------------------------
# Angle wrapping
# ---------------------------------------------------------------------------


def test_wrap_to_pi_range():
    for x in [-3 * np.pi, -np.pi / 2, 0.0, np.pi + 0.1, 5.0, 10.0]:
        w = float(pnt.wrap_to_pi(x))
        assert -np.pi - 1e-9 <= w <= np.pi + 1e-9
        # Wrapping changes the value only by a multiple of 2*pi.
        assert ((x - w) / (2 * np.pi)) == pytest.approx(round((x - w) / (2 * np.pi)), abs=1e-9)


def test_wrap_to_two_pi_range():
    for x in [-1.0, 0.0, 2 * np.pi + 0.3, 10.0]:
        w = float(pnt.wrap_to_two_pi(x))
        assert -1e-9 <= w <= 2 * np.pi + 1e-9


# ---------------------------------------------------------------------------
# Decibel conversions
# ---------------------------------------------------------------------------


def test_decibel_round_trip():
    for value in [0.1, 1.0, 2.0, 100.0]:
        db = pnt.decimal_to_decibel(value)
        assert float(pnt.decibel_to_decimal(db)) == pytest.approx(value, rel=1e-9)


def test_decibel_known_values():
    assert float(pnt.decimal_to_decibel(1.0)) == pytest.approx(0.0, abs=1e-9)
    assert float(pnt.decimal_to_decibel(10.0)) == pytest.approx(10.0, abs=1e-9)
    assert float(pnt.decibel_to_decimal(20.0)) == pytest.approx(100.0, rel=1e-9)


# ---------------------------------------------------------------------------
# Degree / minute / second conversions
# ---------------------------------------------------------------------------


def test_deg_min_sec_round_trip():
    for deg in [0.0, 45.5, 123.456, -30.25]:
        dms = pnt.degrees_to_deg_min_sec(deg)
        back = float(pnt.deg_min_sec_to_degrees(dms))
        assert back == pytest.approx(deg, abs=1e-6)
