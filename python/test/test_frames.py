"""Unit tests for the `pylupnt` reference-frame and ephemeris bindings:
`convert_frame`, `get_frame_center`, `get_body_pos`, `get_body_pos_vel`.

These need the SPICE kernels / EOP data shipped under `LuPNT_data`, so each
test skips cleanly if that data is unavailable.
"""

import numpy as np
import pylupnt as pnt
import pytest

# A representative Earth-orbit epoch and state (LEO-ish, equatorial).
T_TDB = pnt.gregorian_to_time(2026, 1, 1, 0, 0, 0)
RV_ECI = np.array([7000e3, 1000e3, 200e3, 0.5e3, 7.4e3, 0.3e3])
# A Moon-centred state for the lunar frames.
RV_MCI = np.array([3000e3, 0.0, 1500e3, 0.0, 1.4e3, 0.0])


def _skip_if_no_spice():
    try:
        pnt.get_body_pos(T_TDB, pnt.BodyId.EARTH, pnt.BodyId.MOON, pnt.Frame.ECI)
    except Exception as exc:  # pragma: no cover - env dependent
        pytest.skip(f"SPICE/EOP data unavailable: {exc}")


# ---------------------------------------------------------------------------
# convert_frame
# ---------------------------------------------------------------------------


def test_convert_frame_eci_ecef_round_trip():
    _skip_if_no_spice()
    ecef = np.asarray(pnt.convert_frame(T_TDB, RV_ECI, pnt.Frame.ECI, pnt.Frame.ECEF))
    back = np.asarray(pnt.convert_frame(T_TDB, ecef, pnt.Frame.ECEF, pnt.Frame.ECI))
    np.testing.assert_allclose(back, RV_ECI, atol=1e-4)
    # ECEF differs from ECI (Earth has rotated) but conserves geocentric radius.
    assert np.linalg.norm(ecef[:3]) == pytest.approx(np.linalg.norm(RV_ECI[:3]), rel=1e-9)


def test_convert_frame_identity():
    _skip_if_no_spice()
    same = np.asarray(pnt.convert_frame(T_TDB, RV_ECI, pnt.Frame.ECI, pnt.Frame.ECI))
    np.testing.assert_allclose(same, RV_ECI, atol=1e-9)


def test_convert_frame_moon_ci_pa_round_trip():
    _skip_if_no_spice()
    pa = np.asarray(pnt.convert_frame(T_TDB, RV_MCI, pnt.Frame.MOON_CI, pnt.Frame.MOON_PA))
    back = np.asarray(pnt.convert_frame(T_TDB, pa, pnt.Frame.MOON_PA, pnt.Frame.MOON_CI))
    np.testing.assert_allclose(back, RV_MCI, atol=1e-4)
    assert np.linalg.norm(pa[:3]) == pytest.approx(np.linalg.norm(RV_MCI[:3]), rel=1e-9)


def test_convert_frame_vectorized():
    _skip_if_no_spice()
    rv = np.tile(RV_ECI, (4, 1))
    t = T_TDB + np.arange(4) * 600.0
    ecef = np.asarray(pnt.convert_frame(t, rv, pnt.Frame.ECI, pnt.Frame.ECEF))
    assert ecef.shape == (4, 6)


# ---------------------------------------------------------------------------
# get_frame_center
# ---------------------------------------------------------------------------


def test_get_frame_center():
    assert pnt.get_frame_center(pnt.Frame.ECI) == pnt.BodyId.EARTH
    assert pnt.get_frame_center(pnt.Frame.ECEF) == pnt.BodyId.EARTH
    assert pnt.get_frame_center(pnt.Frame.MOON_CI) == pnt.BodyId.MOON
    assert pnt.get_frame_center(pnt.Frame.MOON_PA) == pnt.BodyId.MOON


# ---------------------------------------------------------------------------
# get_body_pos / get_body_pos_vel
# ---------------------------------------------------------------------------


def test_earth_moon_distance():
    _skip_if_no_spice()
    r = np.asarray(pnt.get_body_pos(T_TDB, pnt.BodyId.EARTH, pnt.BodyId.MOON, pnt.Frame.ECI))
    d_km = np.linalg.norm(r) / 1e3
    assert 356_000 < d_km < 407_000  # perigee..apogee of the lunar orbit


def test_earth_sun_distance_is_about_one_au():
    _skip_if_no_spice()
    r = np.asarray(pnt.get_body_pos(T_TDB, pnt.BodyId.EARTH, pnt.BodyId.SUN, pnt.Frame.ECI))
    au_km = np.linalg.norm(r) / 1e3
    assert 1.44e8 < au_km < 1.53e8  # ~1 AU with seasonal variation


def test_get_body_pos_vel_shape_and_consistency():
    _skip_if_no_spice()
    rv = np.asarray(pnt.get_body_pos_vel(T_TDB, pnt.BodyId.EARTH, pnt.BodyId.MOON, pnt.Frame.ECI))
    assert rv.shape == (6,)
    r = np.asarray(pnt.get_body_pos(T_TDB, pnt.BodyId.EARTH, pnt.BodyId.MOON, pnt.Frame.ECI))
    np.testing.assert_allclose(rv[:3], r, atol=1.0)
    # Lunar orbital speed ~1 km/s.
    assert 0.8e3 < np.linalg.norm(rv[3:]) < 1.2e3


def test_get_body_pos_vel_vectorized():
    _skip_if_no_spice()
    t = T_TDB + np.arange(3) * 3600.0
    rv = np.asarray(pnt.get_body_pos_vel(t, pnt.BodyId.EARTH, pnt.BodyId.MOON, pnt.Frame.ECI))
    assert rv.shape == (3, 6)
