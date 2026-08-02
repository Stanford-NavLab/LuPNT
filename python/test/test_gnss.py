"""Unit tests for the `pylupnt` GNSS bindings (`python/bindings/py_gnss.cc`):
`GnssAttitude`, `GnssYawSteering`, `Sp3Loader`, `AntexLoader`,
`RinexNavLoader`, and the shared `GnssConst`/`GnssFreq` enums.

These mirror `cpp/test/agents/test_gnss_attitude.cc`,
`cpp/test/agents/test_gnss_yaw_steering.cc`, and
`cpp/test/interfaces/test_{sp3,antex,rinex_nav}_loader.cc`, exercising the
same C++ classes through their Python bindings -- at minimum, every bound
method is called at least once.
"""

import os

import numpy as np
import pylupnt as pnt
import pytest

ABS_TOL = 1e-6


def _gnss_files_dir():
    """Directory containing the downloaded SP3/BRDC GNSS test files.

    `pnt.get_data_path()` resolves to `<project_root>/data/LuPNT_data`; the
    real GNSS files (downloaded by `scripts/download_gnss_test_fixtures.py` via
    `pnt.Sp3Loader`/`pnt.RinexNavLoader`) live under
    `<project_root>/output/gnss_files`, i.e. two directories up from the data
    path (see `cpp/test/interfaces/test_sp3_loader.cc` for the C++ equivalent
    of this derivation).
    """
    project_root = os.path.dirname(os.path.dirname(pnt.get_data_path()))
    return os.path.join(project_root, "output", "gnss_files")


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


def test_gnss_enums():
    assert pnt.GnssConst.GPS != pnt.GnssConst.GALILEO
    for name in ("GPS", "GLONASS", "GALILEO", "BEIDOU", "QZSS"):
        assert hasattr(pnt.GnssConst, name)

    for name in ("L1", "L2", "L5", "E1", "E6", "E5", "E5a", "E5b"):
        assert hasattr(pnt.GnssFreq, name)
    assert pnt.GnssFreq.L1 != pnt.GnssFreq.L2


# ---------------------------------------------------------------------------
# Antenna
# ---------------------------------------------------------------------------


def test_antenna_gain_pattern():
    try:
        antenna = pnt.Antenna("Block-IIF_ACE")
    except Exception as exc:  # pragma: no cover - env dependent
        pytest.skip(f"Antenna pattern data unavailable: {exc}")

    # Boresight (theta=0) gain is a finite dB value and near the pattern peak.
    g0 = float(antenna.compute_gain(0.0, 0.0))
    assert np.isfinite(g0)

    theta = np.linspace(0.0, np.radians(20.0), 15)
    gains = np.asarray(antenna.compute_gain(theta, 0.0))
    assert gains.shape == theta.shape
    assert np.all(np.isfinite(gains))
    # Peak gain occurs at/near boresight for a nadir-pointing GNSS antenna.
    assert g0 >= np.max(gains) - 1e-6

    gain_matrix = np.asarray(antenna.get_gain_matrix())
    assert gain_matrix.ndim == 2
    n_theta = np.asarray(antenna.get_theta_vector()).shape[0]
    n_phi = np.asarray(antenna.get_phi_vector()).shape[0]
    # The gain grid is indexed by the (theta, phi) sample vectors, in either order.
    assert set(gain_matrix.shape) == {n_theta, n_phi}


# ---------------------------------------------------------------------------
# GnssConstellation
# ---------------------------------------------------------------------------


def test_gnss_constellation_basic():
    const = pnt.GnssConstellation(pnt.GnssConst.GPS)
    assert const.get_gnss_const() == pnt.GnssConst.GPS
    # Freshly constructed: no satellite states loaded yet.
    assert const.get_num_satellites() == 0
    assert const.get_prns() == []


def test_gnss_constellation_set_states():
    const = pnt.GnssConstellation(pnt.GnssConst.GPS)
    t_tai = pnt.convert_time(pnt.gregorian_to_time(2026, 1, 1, 0, 0, 0), pnt.Time.TDB, pnt.Time.TAI)
    t_grid = t_tai + np.arange(3) * 300.0
    prns = [1, 2]
    # Two satellites on trivial straight-line ECI grids [N x 6].
    states = [
        np.tile(np.array([26.56e6, 0, 0, 0, 3.9e3, 0], dtype=float), (3, 1)),
        np.tile(np.array([0, 26.56e6, 0, -3.9e3, 0, 0], dtype=float), (3, 1)),
    ]
    const.set_satellite_states(prns, t_grid, states)
    assert const.get_num_satellites() == 2
    assert sorted(const.get_prns()) == prns
    rv = np.asarray(const.get_satellite_state_eci(1, t_grid[0]))
    assert rv.shape[0] == 6


# ---------------------------------------------------------------------------
# GnssAttitude
# ---------------------------------------------------------------------------


def test_gnss_attitude():
    r_sat = np.array([26560e3, 0.0, 0.0])
    r_sun = np.array([0.3 * pnt.AU, 0.9 * pnt.AU, 0.1 * pnt.AU])

    # Static `compute`
    ex_ref, ey_ref, ez_ref = pnt.GnssAttitude.compute(r_sat, r_sun)

    # Each axis is a unit vector and the triad is orthonormal / right-handed
    for e in (ex_ref, ey_ref, ez_ref):
        assert np.linalg.norm(e) == pytest.approx(1.0, abs=ABS_TOL)
    assert np.dot(ex_ref, ey_ref) == pytest.approx(0.0, abs=ABS_TOL)
    assert np.dot(ey_ref, ez_ref) == pytest.approx(0.0, abs=ABS_TOL)
    np.testing.assert_allclose(ex_ref, np.cross(ey_ref, ez_ref), atol=ABS_TOL)

    # ez points to nadir (along -r_sat)
    np.testing.assert_allclose(ez_ref, -r_sat / np.linalg.norm(r_sat), atol=ABS_TOL)

    # Constructor immediately computes the same triad
    att = pnt.GnssAttitude(r_sat, r_sun)
    np.testing.assert_allclose(att.get_ex(), ex_ref, atol=ABS_TOL)
    np.testing.assert_allclose(att.get_ey(), ey_ref, atol=ABS_TOL)
    np.testing.assert_allclose(att.get_ez(), ez_ref, atol=ABS_TOL)

    # Default-constructed instance starts from the canonical axes
    att_default = pnt.GnssAttitude()
    np.testing.assert_allclose(att_default.get_ex(), [1.0, 0.0, 0.0], atol=ABS_TOL)
    np.testing.assert_allclose(att_default.get_ey(), [0.0, 1.0, 0.0], atol=ABS_TOL)
    np.testing.assert_allclose(att_default.get_ez(), [0.0, 0.0, 1.0], atol=ABS_TOL)

    # `update` (re-)computes and overwrites the cached triad
    att_default.update(r_sat, r_sun)
    np.testing.assert_allclose(att_default.get_ez(), ez_ref, atol=ABS_TOL)

    # Rotation matrix: columns are (ex, ey, ez); orthonormal, det = +1
    R = att.get_rotation_matrix()
    np.testing.assert_allclose(R[:, 0], ex_ref, atol=ABS_TOL)
    np.testing.assert_allclose(R[:, 1], ey_ref, atol=ABS_TOL)
    np.testing.assert_allclose(R[:, 2], ez_ref, atol=ABS_TOL)
    np.testing.assert_allclose(R.T @ R, np.eye(3), atol=ABS_TOL)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=ABS_TOL)

    # Off-boresight angles (theta, phi)
    theta, phi = att.get_angles(att.get_ez())
    assert phi == pytest.approx(0.0, abs=ABS_TOL)

    theta, phi = att.get_angles(att.get_ex())
    assert phi == pytest.approx(pnt.PI_OVER_TWO, abs=ABS_TOL)
    assert theta == pytest.approx(0.0, abs=ABS_TOL)

    theta, phi = att.get_angles(-att.get_ez())
    assert phi == pytest.approx(pnt.PI, abs=ABS_TOL)


def test_gnss_attitude_yaw_steering():
    """Mirrors `agents.gnss_attitude_yaw_steering` in `test_gnss_attitude.cc`:
    verifies that `GnssAttitude` is connected to the implemented
    `GnssYawSteering` nominal yaw-steering law (Eq. 1, Cheng et al., 2025),
    https://doi.org/10.1016/j.asr.2024.10.064 -- i.e. the new
    velocity-aware `compute`/`compute_from_yaw_angle`/`update` overloads
    reproduce exactly the existing Sun-direction-based attitude frame.
    """
    ys = pnt.GnssYawSteering

    r_sat = np.array([26560e3, 0.0, 0.0])
    v_sat = np.array([0.0, 3873.7, 0.0])
    r_sun = np.array([0.3 * pnt.AU, 0.9 * pnt.AU, 0.1 * pnt.AU])

    # Looser tolerance for cross-checking two different floating-point paths
    # (geometric Sun-pointing vs. trig-chain through GnssYawSteering) that are
    # mathematically identical but numerically distinct.
    cross_path_tol = 1e-7

    ex_geom, ey_geom, ez_geom = pnt.GnssAttitude.compute(r_sat, r_sun)

    # `compute(r_sat, v_sat, r_sun)` routes through GnssYawSteering's nominal
    # yaw-steering law and matches the Sun-direction-based `compute`
    ex_yaw, ey_yaw, ez_yaw = pnt.GnssAttitude.compute(r_sat, v_sat, r_sun)
    np.testing.assert_allclose(ex_yaw, ex_geom, atol=cross_path_tol)
    np.testing.assert_allclose(ey_yaw, ey_geom, atol=cross_path_tol)
    np.testing.assert_allclose(ez_yaw, ez_geom, atol=cross_path_tol)

    # `compute_from_yaw_angle` with the nominal yaw angle phi_nom =
    # nominal_yaw_angle(beta, mu) reproduces the Sun-pointing frame exactly
    beta = ys.beta_angle(r_sat, v_sat, r_sun)
    mu = ys.orbit_angle(r_sat, v_sat, r_sun)
    phi_nom = ys.nominal_yaw_angle(beta, mu)

    ex_phi, ey_phi, ez_phi = pnt.GnssAttitude.compute_from_yaw_angle(r_sat, v_sat, phi_nom)
    np.testing.assert_allclose(ex_phi, ex_geom, atol=cross_path_tol)
    np.testing.assert_allclose(ey_phi, ey_geom, atol=cross_path_tol)
    np.testing.assert_allclose(ez_phi, ez_geom, atol=cross_path_tol)

    # `compute_from_yaw_angle` produces an orthonormal, right-handed,
    # nadir-pointing triad for an arbitrary (non-nominal) yaw angle too
    ex, ey, ez = pnt.GnssAttitude.compute_from_yaw_angle(r_sat, v_sat, 0.37)
    for e in (ex, ey, ez):
        assert np.linalg.norm(e) == pytest.approx(1.0, abs=ABS_TOL)
    assert np.dot(ex, ey) == pytest.approx(0.0, abs=ABS_TOL)
    assert np.dot(ey, ez) == pytest.approx(0.0, abs=ABS_TOL)
    np.testing.assert_allclose(ex, np.cross(ey, ez), atol=ABS_TOL)
    np.testing.assert_allclose(ez, -r_sat / np.linalg.norm(r_sat), atol=ABS_TOL)

    # `update(r_sat, v_sat, r_sun)` caches the same triad as the
    # yaw-law-driven `compute`
    att = pnt.GnssAttitude()
    att.update(r_sat, v_sat, r_sun)
    np.testing.assert_allclose(att.get_ex(), ex_yaw, atol=ABS_TOL)
    np.testing.assert_allclose(att.get_ey(), ey_yaw, atol=ABS_TOL)
    np.testing.assert_allclose(att.get_ez(), ez_yaw, atol=ABS_TOL)


# ---------------------------------------------------------------------------
# GnssYawSteering
# ---------------------------------------------------------------------------
#
# Mirrors `cpp/test/agents/test_gnss_yaw_steering.cc`: exercises the GPS/
# Galileo/BDS-3 yaw-attitude steering laws of Cheng et al. (2025),
# https://doi.org/10.1016/j.asr.2024.10.064 (see `GnssYawSteering` in
# `lupnt/agents/gnss_yaw_steering.h` for the equation each method mirrors).


def test_gnss_yaw_steering_geometry():
    ys = pnt.GnssYawSteering
    r = 26560e3
    v = np.sqrt(pnt.GM_EARTH / r)
    r_sat = np.array([r, 0.0, 0.0])
    v_sat = np.array([0.0, v, 0.0])

    # beta: Sun in the orbital plane -> 0; above/below -> +/- 90 deg
    assert ys.beta_angle(
        r_sat, v_sat, np.array([0.3 * pnt.AU, 0.4 * pnt.AU, 0.0])
    ) == pytest.approx(0.0, abs=ABS_TOL)
    assert ys.beta_angle(r_sat, v_sat, np.array([0.0, 0.0, pnt.AU])) == pytest.approx(
        pnt.PI_OVER_TWO, abs=ABS_TOL
    )
    assert ys.beta_angle(r_sat, v_sat, np.array([0.0, 0.0, -pnt.AU])) == pytest.approx(
        -pnt.PI_OVER_TWO, abs=ABS_TOL
    )

    # mu: zero at the midnight point, +/- pi at the noon point
    assert ys.orbit_angle(r_sat, v_sat, np.array([-pnt.AU, 0.0, 0.0])) == pytest.approx(
        0.0, abs=ABS_TOL
    )
    assert abs(ys.orbit_angle(r_sat, v_sat, np.array([pnt.AU, 0.0, 0.0]))) == pytest.approx(
        pnt.PI, abs=ABS_TOL
    )

    # eta = wrap_to_pi(mu - pi)
    assert ys.orbit_noon_angle(0.0) == pytest.approx(-pnt.PI, abs=ABS_TOL)
    assert ys.orbit_noon_angle(pnt.PI) == pytest.approx(0.0, abs=ABS_TOL)


def test_gnss_yaw_steering_nominal_model():
    ys = pnt.GnssYawSteering
    mu = np.deg2rad(30.0)

    # Eq. (1): sign(phi) is always opposite sign(beta); symmetric in beta -> -beta
    phi_pos = ys.nominal_yaw_angle(np.deg2rad(5.0), mu)
    phi_neg = ys.nominal_yaw_angle(np.deg2rad(-5.0), mu)
    assert phi_pos < 0.0
    assert phi_neg > 0.0
    assert phi_pos == pytest.approx(-phi_neg, abs=ABS_TOL)

    # At mu = 0 / pi (midnight / noon points), |phi| = 90 deg
    assert abs(ys.nominal_yaw_angle(np.deg2rad(5.0), 0.0)) == pytest.approx(
        pnt.PI_OVER_TWO, abs=ABS_TOL
    )

    # Eq. (2): nominal yaw rate matches a numerical derivative of Eq. (1) w.r.t. time
    mu_dot = pnt.TWO_PI / 43200.0
    dmu = 1e-6
    dphi_dmu = (
        ys.nominal_yaw_angle(np.deg2rad(5.0), mu + dmu)
        - ys.nominal_yaw_angle(np.deg2rad(5.0), mu - dmu)
    ) / (2 * dmu)
    assert ys.nominal_yaw_rate(np.deg2rad(5.0), mu, mu_dot) == pytest.approx(
        dphi_dmu * mu_dot, abs=1e-6
    )


def test_gnss_yaw_steering_gps():
    ys = pnt.GnssYawSteering
    ts, te = 1000.0, 4000.0
    phi_ts = np.deg2rad(170.0)
    phi_te = np.deg2rad(-170.0)

    # GPS Block IIF shadow crossing: linear interpolation  [Eq. 3-4]
    assert ys.gps_iif_shadow_yaw_angle(ts, ts, te, phi_ts, phi_te) == pytest.approx(
        phi_ts, abs=ABS_TOL
    )
    assert ys.gps_iif_shadow_yaw_angle(te, ts, te, phi_ts, phi_te) == pytest.approx(
        phi_te, abs=ABS_TOL
    )
    t_mid = 0.5 * (ts + te)
    assert ys.gps_iif_shadow_yaw_angle(t_mid, ts, te, phi_ts, phi_te) == pytest.approx(
        0.5 * (phi_ts + phi_te), abs=ABS_TOL
    )

    # GPS Block IIF noon turn: maximum hardware yaw rate R = 0.11 deg/s  [Eq. 5]
    R = np.deg2rad(0.11)
    t = ts + 100.0
    phi = ys.gps_iif_noon_turn_yaw_angle(t, ts, phi_ts, np.deg2rad(1.0))
    assert phi == pytest.approx(phi_ts - R * (t - ts), abs=ABS_TOL)
    assert phi < phi_ts

    # GPS Block IIR midnight/noon turns: hardware yaw rate R = 0.20 deg/s  [Eq. 6-7]
    R = np.deg2rad(0.20)
    beta = np.deg2rad(0.3)
    assert ys.gps_iir_midnight_turn_yaw_angle(t, ts, phi_ts, beta) == pytest.approx(
        phi_ts + R * (t - ts), abs=ABS_TOL
    )
    assert ys.gps_iir_noon_turn_yaw_angle(t, ts, phi_ts, beta) == pytest.approx(
        phi_ts - R * (t - ts), abs=ABS_TOL
    )


def test_gnss_yaw_steering_galileo():
    ys = pnt.GnssYawSteering

    # Sun reference vector is a unit vector matching Eq. (9) component-wise
    eta, beta = np.deg2rad(20.0), np.deg2rad(3.0)
    S = ys.galileo_sun_vector(eta, beta)
    np.testing.assert_allclose(np.linalg.norm(S), 1.0, atol=1e-12)
    np.testing.assert_allclose(
        S, [np.sin(eta) * np.cos(beta), np.sin(beta), np.cos(eta) * np.cos(beta)], atol=ABS_TOL
    )

    # IOV nominal model (Eq. 8-9) is *algebraically identical* to the generic
    # nominal model (Eq. 1) once eta = orbit_noon_angle(mu) = mu - pi is
    # substituted: -S_y = -sin(beta), -S_x = -sin(eta)*cos(beta)
    # = sin(mu)*cos(beta), so (dropping the common positive scale factor)
    # ATAN2(-S_y, -S_x) = ATAN2(-sin(beta), sin(mu)*cos(beta))
    #                   = ATAN2(-tan(beta), sin(mu)) = Eq. (1).
    mu = np.deg2rad(50.0)
    eta = ys.orbit_noon_angle(mu)
    assert ys.galileo_iov_nominal_yaw_angle(eta, beta) == pytest.approx(
        ys.nominal_yaw_angle(beta, mu), abs=1e-9
    )

    # IOV eclipse model (Eq. 10-11) returns a finite, valid yaw angle
    phi_nom = ys.nominal_yaw_angle(beta, mu)
    phi_eclipse = ys.galileo_iov_eclipse_yaw_angle(eta, beta, phi_nom)
    assert np.isfinite(phi_eclipse)
    assert -pnt.PI <= phi_eclipse <= pnt.PI

    # FOC yaw-steering law (Eq. 12): half-cosine oscillation between phi_s and
    # 90 deg * SIGN(1, phi_s), with period 5656 s
    ts = 0.0
    phi_s = np.deg2rad(150.0)
    period = 5656.0
    assert ys.galileo_foc_yaw_angle(ts, ts, phi_s) == pytest.approx(phi_s, abs=ABS_TOL)
    assert ys.galileo_foc_yaw_angle(ts + period, ts, phi_s) == pytest.approx(phi_s, abs=ABS_TOL)
    expected_half = np.pi * ys.sign(1.0, phi_s) - phi_s
    assert ys.galileo_foc_yaw_angle(ts + period / 2.0, ts, phi_s) == pytest.approx(
        expected_half, abs=ABS_TOL
    )


def test_gnss_yaw_steering_gps3():
    ys = pnt.GnssYawSteering

    # Sun direction unit vector is a unit vector matching Eq. (11) component-wise
    beta, mu = np.deg2rad(3.0), np.deg2rad(20.0)
    s = ys.gps3_sun_vector(beta, mu)
    np.testing.assert_allclose(np.linalg.norm(s), 1.0, atol=1e-12)
    np.testing.assert_allclose(
        s, [np.cos(beta) * np.sin(mu), -np.sin(beta), np.cos(beta) * np.cos(mu)], atol=ABS_TOL
    )

    # The new (mu-direct) Sun-vector parameterization yields the *same* nominal
    # yaw angle as the generic nominal model (Eq. 12 reduces algebraically to
    # Eq. 1, just as the eta-based Galileo parameterization does):
    #   ATAN2(s_y, s_x) = ATAN2(-sin(beta), cos(beta)*sin(mu))
    #                   = ATAN2(-tan(beta), sin(mu)) = Eq. (1)
    beta2, mu2 = np.deg2rad(1.0), np.deg2rad(50.0)
    s2 = ys.gps3_sun_vector(beta2, mu2)
    assert np.arctan2(s2[1], s2[0]) == pytest.approx(ys.nominal_yaw_angle(beta2, mu2), abs=1e-9)

    # Improved eclipse model (Eq. 13, 17-18) returns a finite, valid yaw angle
    phi_nom = ys.nominal_yaw_angle(beta2, mu2)
    phi_eclipse = ys.gps3_eclipse_yaw_angle(beta2, mu2, phi_nom)
    assert np.isfinite(phi_eclipse)
    assert -pnt.PI <= phi_eclipse <= pnt.PI

    # At the orbit noon/midnight point (s_x = 0), the cosine weighting g = 1,
    # so Eq. (18) collapses to the full collinearity-region clamp
    # s_y* = SIGN(1, phi_nom) * gamma_y  (gamma_y = sin(5.8 deg) for GPS III)
    mu0 = 0.0
    phi_nom0 = ys.nominal_yaw_angle(beta2, mu0)
    phi0 = ys.gps3_eclipse_yaw_angle(beta2, mu0, phi_nom0)
    gamma_y = np.sin(np.deg2rad(5.8))
    expected0 = np.arctan2(ys.sign(1.0, phi_nom0) * gamma_y, 0.0)
    assert phi0 == pytest.approx(expected0, abs=1e-9)


def test_gnss_yaw_steering_bds3():
    ys = pnt.GnssYawSteering

    # CAST WHU model (Eq. 13-14): same form as Eq. 12 with different periods
    ts = 200.0
    phi_s = np.deg2rad(-160.0)
    assert ys.bds3_cast_igso_yaw_angle(ts, ts, phi_s) == pytest.approx(phi_s, abs=ABS_TOL)
    assert ys.bds3_cast_meo_yaw_angle(ts, ts, phi_s) == pytest.approx(phi_s, abs=ABS_TOL)
    assert ys.bds3_cast_igso_yaw_angle(ts + 5740.0, ts, phi_s) == pytest.approx(phi_s, abs=ABS_TOL)
    assert ys.bds3_cast_meo_yaw_angle(ts + 3090.0, ts, phi_s) == pytest.approx(phi_s, abs=ABS_TOL)
    # The two models differ at intermediate epochs (different oscillation periods)
    t = ts + 1000.0
    assert ys.bds3_cast_igso_yaw_angle(t, ts, phi_s) != pytest.approx(
        ys.bds3_cast_meo_yaw_angle(t, ts, phi_s), abs=1e-3
    )

    # SECM CSNO model (Eq. 15): the yaw target sign flips discontinuously across beta = 0
    mu = np.deg2rad(10.0)
    three_deg = np.deg2rad(3.0)
    phi_pos = ys.bds3_secm_csno_yaw_angle(np.deg2rad(1.0), mu)
    phi_neg = ys.bds3_secm_csno_yaw_angle(np.deg2rad(-1.0), mu)
    assert phi_pos == pytest.approx(np.arctan2(-np.tan(three_deg), np.sin(mu)), abs=ABS_TOL)
    assert phi_neg == pytest.approx(np.arctan2(np.tan(three_deg), np.sin(mu)), abs=ABS_TOL)
    assert abs(phi_pos - phi_neg) > 1e-3

    # SECM MCSNO model (Eq. 16): matches CSNO before t0 (no discontinuity test
    # here -- that's the point of the model -- just check finiteness/continuity)
    t0, ts2, te = 1000.0, 1200.0, 1380.0
    mu_dot = pnt.TWO_PI / 43200.0
    mu_ts = np.deg2rad(90.0)
    beta_dot = -2e-7

    phi_before = ys.bds3_secm_mcsno_yaw_angle(
        500.0, t0, ts2, te, np.deg2rad(0.5), 1e-7, np.deg2rad(100.0), 0.0, mu_dot, 0.0
    )
    assert phi_before == pytest.approx(
        np.arctan2(-np.tan(three_deg), np.sin(np.deg2rad(100.0))), abs=ABS_TOL
    )

    # Per Eq. (16), `phi_ts` is "the nominal yaw angle at the [transition]
    # starting time ts" -- i.e. it is supplied by the caller (typically the
    # [t0, ts) branch's value evaluated at mu = mu_ts), not computed
    # internally. At t = ts2 the function lands in the [ts, te) branch (the
    # `t < ts` test is false at equality) and the linear-transition term
    # `(mu - mu_ts)` vanishes, so it returns the supplied `phi_ts` unchanged --
    # we compute the continuity-enforcing boundary value directly and supply
    # it, mirroring how an upstream caller would chain the two branches.
    expected_phi_ts = np.arctan2(np.tan(ys.sign(three_deg, beta_dot)), np.sin(mu_ts))
    phi_at_ts = ys.bds3_secm_mcsno_yaw_angle(
        ts2, t0, ts2, te, 0.0, beta_dot, mu_ts, mu_ts, mu_dot, expected_phi_ts
    )
    assert phi_at_ts == pytest.approx(expected_phi_ts, abs=ABS_TOL)

    # Linear-transition rate matches the documented 0.055 deg/s
    transition_rate = np.deg2rad(0.055)
    mu_t = mu_ts + mu_dot * 60.0
    phi_t = ys.bds3_secm_mcsno_yaw_angle(
        ts2 + 60.0, t0, ts2, te, 0.0, beta_dot, mu_t, mu_ts, mu_dot, expected_phi_ts
    )
    assert phi_t == pytest.approx(
        expected_phi_ts - ys.sign(transition_rate, beta_dot) * 60.0, abs=1e-9
    )


def test_gnss_yaw_steering_sign():
    ys = pnt.GnssYawSteering
    # FORTRAN SIGN(a, b): magnitude of `a` with the sign of `b`
    assert ys.sign(2.0, 3.0) == pytest.approx(2.0, abs=ABS_TOL)
    assert ys.sign(2.0, -3.0) == pytest.approx(-2.0, abs=ABS_TOL)
    assert ys.sign(-2.0, 3.0) == pytest.approx(2.0, abs=ABS_TOL)
    assert ys.sign(-2.0, -3.0) == pytest.approx(-2.0, abs=ABS_TOL)
    assert ys.sign(-2.0, 0.0) == pytest.approx(2.0, abs=ABS_TOL)


# ---------------------------------------------------------------------------
# Sp3Loader
# ---------------------------------------------------------------------------


def test_sp3_loader():
    sp3_dir = os.path.join(_gnss_files_dir(), "sp3")
    sp3_file_1 = os.path.join(sp3_dir, "COD0MGXFIN_20260140000_01D_05M_ORB.SP3")
    sp3_file_2 = os.path.join(sp3_dir, "COD0MGXFIN_20260150000_01D_05M_ORB.SP3")
    sat_id = "G01"

    if not os.path.isfile(sp3_file_1):
        pytest.skip(f"SP3 test data not found at {sp3_file_1}")

    # Default constructor: nothing loaded
    loader_empty = pnt.Sp3Loader()
    assert loader_empty.get_satellites() == []
    assert not loader_empty.has_satellite(sat_id)

    # Single-file constructor
    loader = pnt.Sp3Loader(sp3_file_1)
    assert len(loader.get_satellites()) > 0
    assert loader.has_satellite(sat_id)
    assert not loader.has_satellite("X99")

    t_min, t_max = loader.get_time_span(sat_id)
    assert t_max > t_min
    t_mid = 0.5 * (t_min + t_max)

    rv_ecef, clock_bias_s = loader.get_pos_vel_clock(sat_id, t_mid)
    r = np.linalg.norm(rv_ecef[:3])
    v = np.linalg.norm(rv_ecef[3:])
    assert 2.0e7 < r < 3.0e7  # GPS geocentric radius ~26,560 km
    assert 1.0e3 < v < 1.0e4  # GPS orbital speed ~3.9 km/s
    assert abs(clock_bias_s) < 1e-3

    rv_only = loader.get_pos_vel(sat_id, t_mid)
    np.testing.assert_allclose(rv_only, rv_ecef, atol=1e-6)

    # `load_file` / multi-file constructor merge & extend the time span
    loader.load_file(sp3_file_2)
    t2_min, t2_max = loader.get_time_span(sat_id)
    assert t2_min <= t_min
    assert t2_max > t_max

    loader_multi = pnt.Sp3Loader([sp3_file_1, sp3_file_2])
    t3_min, t3_max = loader_multi.get_time_span(sat_id)
    assert t3_min == pytest.approx(t2_min, abs=1e-6)
    assert t3_max == pytest.approx(t2_max, abs=1e-6)


# ---------------------------------------------------------------------------
# SP3 / RINEX CDDIS download helpers (filename / URL builders)
# ---------------------------------------------------------------------------


def _tai(y, mo, d, h=12, mi=0, s=0.0):
    t = pnt.gregorian_to_time(y, mo, d, h, mi, s)
    return pnt.convert_time(t, pnt.Time.TDB, pnt.Time.TAI)


def test_sp3_download_helpers_filename_and_url():
    """Sp3Loader static builders map an epoch to the COD MGEX product name/URL."""
    t = _tai(2026, 1, 14)
    fname = pnt.Sp3Loader.filename_for_epoch(t, pnt.Time.TAI)
    assert fname == "COD0MGXFIN_20260140000_01D_05M_ORB.SP3"

    url = pnt.Sp3Loader.url_for_epoch(t, pnt.Time.TAI)
    assert url.startswith("https://cddis.nasa.gov/archive/gnss/products/")
    assert url.endswith("COD0MGXFIN_20260140000_01D_05M_ORB.SP3.gz")

    # The default time scale is UTC; a UTC-built epoch resolves to the same day.
    t_utc = pnt.gregorian_to_time(2026, 1, 15, 0, 0, 0)
    assert pnt.Sp3Loader.filename_for_epoch(t_utc) == "COD0MGXFIN_20260150000_01D_05M_ORB.SP3"


def test_rinex_download_helper_filename():
    """RinexNavLoader static builder maps an epoch to the BRDC product name."""
    t = _tai(2026, 1, 14)
    fname = pnt.RinexNavLoader.filename_for_epoch(t, pnt.Time.TAI)
    assert fname == "BRDC00IGS_R_20260140000_01D_MN.rnx"


def test_download_file_for_epoch_returns_cached_path():
    """When the product is already cached, download_file_for_epoch returns that
    path without hitting the network (so it runs offline once fixtures exist)."""
    sp3_file = os.path.join(_gnss_files_dir(), "sp3", "COD0MGXFIN_20260140000_01D_05M_ORB.SP3")
    if not os.path.isfile(sp3_file):
        pytest.skip(f"SP3 fixture not cached at {sp3_file}")

    t = _tai(2026, 1, 14)
    path = pnt.Sp3Loader.download_file_for_epoch(t, pnt.Time.TAI)
    assert os.path.isfile(path)
    assert os.path.basename(path) == "COD0MGXFIN_20260140000_01D_05M_ORB.SP3"

    brdc_file = os.path.join(_gnss_files_dir(), "brdc", "BRDC00IGS_R_20260140000_01D_MN.rnx")
    if os.path.isfile(brdc_file):
        bpath = pnt.RinexNavLoader.download_file_for_epoch(t, pnt.Time.TAI)
        assert os.path.isfile(bpath)
        assert os.path.basename(bpath) == "BRDC00IGS_R_20260140000_01D_MN.rnx"


# ---------------------------------------------------------------------------
# AntexLoader
# ---------------------------------------------------------------------------


def test_antex_loader():
    antex_file = os.path.join(pnt.get_data_path(), "gnss", "igs20.atx")
    if not os.path.isfile(antex_file):
        pytest.skip(f"ANTEX test data not found at {antex_file}")

    t_tai = pnt.gregorian_to_time(2026, 1, 14, 0, 0, 0)
    # SP3 (center-of-mass) ECEF position sample for G01 at the same epoch.
    r_sat_ecef = np.array([21175.826427e3, 13096.416599e3, 9341.310042e3])

    # Default constructor: nothing loaded
    antex_empty = pnt.AntexLoader()
    assert not antex_empty.has_satellite(pnt.GnssConst.GPS, 1)

    # Single-file constructor
    antex = pnt.AntexLoader(antex_file)
    assert antex.has_satellite(pnt.GnssConst.GPS, 1)
    assert not antex.has_satellite(pnt.GnssConst.GPS, 99)

    codes = antex.get_available_freq_codes(pnt.GnssConst.GPS, 1, t_tai)
    assert len(codes) > 0
    assert "G01" in codes

    # Both `get_pco` overloads agree
    pco_l1 = antex.get_pco(pnt.GnssConst.GPS, 1, pnt.GnssFreq.L1, t_tai)
    pco_g01 = antex.get_pco("G", 1, "G01", t_tai)
    np.testing.assert_allclose(pco_l1, pco_g01, atol=1e-6)
    assert 0.1 < np.linalg.norm(pco_l1) < 100.0  # PCO magnitude on the order of meters

    # `load_file` merges additional satellite entries
    antex_loaded = pnt.AntexLoader()
    assert not antex_loaded.has_satellite(pnt.GnssConst.GPS, 1)
    antex_loaded.load_file(antex_file)
    assert antex_loaded.has_satellite(pnt.GnssConst.GPS, 1)

    # Static helpers: IJK rotation, PCO correction, identifiers
    # Note: Cijk is intentionally not orthonormal -- it is a verbatim port of
    # `phase_center_offset.py::ijk_to_ecef_rot` (see AntexLoader::ComputeIjkToEcefRotation).
    Cijk = pnt.AntexLoader.compute_ijk_to_ecef_rotation(t_tai, r_sat_ecef)
    np.testing.assert_allclose(Cijk[:, 2], -r_sat_ecef / np.linalg.norm(r_sat_ecef), atol=1e-6)
    np.testing.assert_allclose(Cijk[:, 0], np.cross(Cijk[:, 1], Cijk[:, 2]), atol=1e-6)

    corrected = pnt.AntexLoader.apply_pco_correction_ecef(t_tai, r_sat_ecef, pco_l1)
    expected = r_sat_ecef + Cijk @ pco_l1
    np.testing.assert_allclose(corrected, expected, atol=1e-6)
    assert np.linalg.norm(corrected - r_sat_ecef) > 0.0

    assert pnt.AntexLoader.gnss_letter(pnt.GnssConst.GPS) == "G"
    assert pnt.AntexLoader.gnss_letter(pnt.GnssConst.GLONASS) == "R"
    assert pnt.AntexLoader.gnss_letter(pnt.GnssConst.GALILEO) == "E"
    assert pnt.AntexLoader.gnss_letter(pnt.GnssConst.BEIDOU) == "C"
    assert pnt.AntexLoader.gnss_letter(pnt.GnssConst.QZSS) == "J"

    assert pnt.AntexLoader.sat_id(pnt.GnssConst.GPS, 1) == "G01"
    assert pnt.AntexLoader.sat_id(pnt.GnssConst.GPS, 32) == "G32"
    assert pnt.AntexLoader.sat_id(pnt.GnssConst.GALILEO, 11) == "E11"
    assert pnt.AntexLoader.sat_id(pnt.GnssConst.BEIDOU, 6) == "C06"


# ---------------------------------------------------------------------------
# RinexNavLoader
# ---------------------------------------------------------------------------


def test_rinex_nav_loader():
    brdc_dir = os.path.join(_gnss_files_dir(), "brdc")
    nav_file = os.path.join(brdc_dir, "BRDC00IGS_R_20260140000_01D_MN.rnx")
    sat_id = "G01"

    if not os.path.isfile(nav_file):
        pytest.skip(f"RINEX nav test data not found at {nav_file}")

    t_tai = pnt.gregorian_to_time(2026, 1, 14, 1, 0, 0)

    # Default constructor: nothing loaded
    loader_empty = pnt.RinexNavLoader()
    assert loader_empty.get_satellites() == []
    assert not loader_empty.has_satellite(sat_id)

    # Single-file constructor
    loader = pnt.RinexNavLoader(nav_file)
    sats = loader.get_satellites()
    assert len(sats) > 0
    assert sat_id in sats
    assert not loader.has_satellite("X99")
    # GLONASS ('R', tabulated state-vector ephemeris) is intentionally excluded
    assert all(not s.startswith("R") for s in sats)

    rv_ecef, clock_corr_s = loader.get_pos_vel_clock(sat_id, t_tai)
    r = np.linalg.norm(rv_ecef[:3])
    v = np.linalg.norm(rv_ecef[3:])
    assert 2.0e7 < r < 3.0e7
    assert 1.0e3 < v < 1.0e4
    assert abs(clock_corr_s) < 1e-3

    rv_only = loader.get_pos_vel(sat_id, t_tai)
    np.testing.assert_allclose(rv_only, rv_ecef, atol=1e-6)

    # `load_file` / multi-file constructor accumulate navigation messages
    loader.load_file(nav_file)
    assert loader.has_satellite(sat_id)

    loader_multi = pnt.RinexNavLoader([nav_file])
    rv_multi = loader_multi.get_pos_vel(sat_id, t_tai)
    np.testing.assert_allclose(rv_multi, rv_only, atol=1e-6)
