import numpy as np
import pylupnt as pnt


def load():
    # Epoch(TAI)
    t0_tai_str = "2022-01-15T00:00:00"
    N_sat, N_planes = 1, 1

    # Classical orbital elements (a, e, i, W, w, M) [km, -, rad, rad, rad, rad]
    sma = 5740  # [km] Semi-major axis
    ecc = [0.58]  # [-] Eccentricity
    inc = np.deg2rad([54.856])  # [rad] Inclination
    raan = np.deg2rad([0])  # [rad] Right ascension of the ascending node
    aop = np.deg2rad([86.322])  # [rad] Argument of periapsis
    ma = np.deg2rad([0])  # [rad] Mean anomaly

    t0_tai = pnt.spice.string2tai(t0_tai_str)
    coe_op = np.zeros((N_sat, 6))
    N_sat_plane = N_sat // N_planes
    for i_pl in range(N_planes):
        for i_spl in range(N_sat_plane):
            i_sat = i_pl * N_sat_plane + i_spl
            coe_op[i_sat] = np.array(
                [sma, ecc[i_pl], inc[i_pl], raan[i_pl], aop[i_pl], ma[i_spl]]
            )

    rv0_m2sc_op = pnt.classical2cart(coe_op, pnt.GM_MOON)
    rv0_m2sc_ci = pnt.convert_frame(t0_tai, rv0_m2sc_op, pnt.MOON_OP, pnt.MOON_CI)

    # Time
    sma = coe_op[0, 0]  # [km] Semi-major axis
    period = 2 * np.pi * np.sqrt(np.power(sma, 3) / pnt.GM_MOON)  # [s] Orbital period
    Dt = 5 * pnt.SECS_MINUTE  # [s] Simulation time step
    dt = 5 * pnt.SECS_MINUTE  # [s] Propagation time step
    tf = period  # [s] Simulation final time
    Nt = int(tf / Dt)  # [-] Number of time steps
    tspan = np.linspace(0, tf, Nt)  # [s] Time since first epoch
    t_tai = t0_tai + tspan  # [s] Epochs (TAI)

    # Dynamics (three-body)
    dyn = pnt.NBodyDynamics()
    dyn.add_body(pnt.Body.Moon())
    dyn.add_body(pnt.Body.Earth())
    dyn.set_frame(pnt.MOON_CI)
    dyn.set_time_step(dt)

    # Propagation
    # rv_from2to_frame [km, km/s] (x, y, z, vx, vy, vz)
    rv_m2sc_ci = np.zeros((N_sat, Nt, 6))
    rv_m2sc_pa = np.zeros((N_sat, Nt, 6))
    for i_sat in range(N_sat):
        rv_m2sc_ci[i_sat] = dyn.propagate(rv0_m2sc_ci[i_sat], t0_tai, t_tai)
        rv_m2sc_pa[i_sat] = pnt.convert_frame(
            t_tai, rv_m2sc_ci[i_sat], pnt.MOON_CI, pnt.MOON_PA
        )

    rv_m2e_ci = pnt.spice.get_body_pos_vel(t_tai, pnt.MOON, pnt.EARTH)
    rv_m2e_pa = pnt.convert_frame(t_tai, rv_m2e_ci, pnt.MOON_CI, pnt.MOON_PA)
    rv_m2s_ci = pnt.spice.get_body_pos_vel(t_tai, pnt.MOON, pnt.SUN)
    rv_m2s_pa = pnt.convert_frame(t_tai, rv_m2s_ci, pnt.MOON_CI, pnt.MOON_PA)

    # Directions
    e_sc2m = np.array(-rv_m2sc_ci[:, :, 0:3])
    e_sc2s = np.array(rv_m2s_ci[None, :, 0:3] - rv_m2sc_ci[:, :, 0:3])
    e_sc2e = np.array(rv_m2e_ci[None, :, 0:3] - rv_m2sc_ci[:, :, 0:3])
    e_sc2m /= np.linalg.norm(e_sc2m, axis=2)[:, :, None]
    e_sc2s /= np.linalg.norm(e_sc2s, axis=2)[:, :, None]
    e_sc2e /= np.linalg.norm(e_sc2e, axis=2)[:, :, None]

    # Yaw-Steering
    ez_sc = e_sc2m
    ey_sc = pnt.cross_norm(e_sc2m, e_sc2s)
    ex_sc = pnt.cross_norm(ey_sc, ez_sc)

    # Sun angle
    sun_angle_cos = np.sum(e_sc2s * ex_sc, axis=2)
    sun_angle = np.arccos(np.clip(sun_angle_cos, -1, 1))

    # Passive rotation matrices
    R_mi2sc = np.zeros((N_sat, Nt, 3, 3))
    R_mi2sc[:, :, 0, :] = ex_sc
    R_mi2sc[:, :, 1, :] = ey_sc
    R_mi2sc[:, :, 2, :] = ez_sc

    assert np.allclose(np.linalg.norm(ex_sc, axis=-1), 1)
    assert np.allclose(np.linalg.norm(ey_sc, axis=-1), 1)
    assert np.allclose(np.linalg.norm(ez_sc, axis=-1), 1)

    R_pa2mi = np.zeros((Nt, 3, 3))
    for i in range(Nt):
        R_pa2mi[i] = pnt.spice.get_frame_conversion_matrix(
            t_tai[i], pnt.MOON_PA, pnt.MOON_CI
        )[:3, :3]

    # Spacecraft frame
    cam_dir = np.array([0, 0, 1])
    cam_up = np.array([1, 0, 0])
    cam_right = np.array([0, 1, 0])
    ez_c = cam_dir
    ey_c = -cam_up
    ex_c = cam_right
    R_sc2ocv = np.array([ex_c, ey_c, ez_c])
