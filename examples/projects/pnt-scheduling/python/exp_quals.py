import os
import util

import numpy as np
import pylupnt as pnt
import pathfinder_data
from math import ceil, floor
from problem import Request, ServiceWindow, PntSchedulingProblem
from solvers import (
    SmdpForwardSearchSolver,
    SmdpMctsSolver,
    DiscreteTimeIpSolver,
    RuleBasedSolver,
)
import logging

from functools import partial

cache_path = "cache/quals/"

solvers = {
    "RB": RuleBasedSolver,
    "FS": SmdpForwardSearchSolver,
    "MCTS": SmdpMctsSolver,
    "IP": DiscreteTimeIpSolver,
}


def solve_wrapper(
    config: dict, problem: PntSchedulingProblem, solver: str, n_jobs: int
) -> list[dict]:
    """
    Load or recompute results for a solver and specific parameters.

    Args:
        config (dict): Configuration including date, duration factor, and solver parameters
        problem (PntSchedulingProblem): Problem
        solver (str): Solver name
        n_jobs (int): Number of jobs

    Returns:
        list[dict]: Results
    """
    # Solve function
    func = partial(
        util.solve, problem=problem, solver=solvers[solver], seed=config["solver_seed"]
    )

    # Iterable with all the parameter combinations
    it = util.iterate_params(config[solver])

    # Filepath for caching
    filepath = os.path.join(
        cache_path,
        solver,
        util.create_hash(config["date"], config["duration_factor"], config[solver]),
    )

    # Load or recompute results
    results = util.load_or_recompute(
        filepath, util.run_parallel, func, it, n_jobs=n_jobs, desc=solver
    )

    # Sort results by duration percentage and computing time
    results.sort(key=lambda x: 1e3 * x["total"] - x["time"], reverse=True)
    return results


def solve_problem(config: dict, n_jobs: int = 5) -> dict:
    """
    Solve the problem using all solvers.

    Args:
        config (dict): Configuration including date, duration factor, and solver parameters
        n_jobs (int): Number of jobs

    Returns:
        dict: Results
    """

    # Create the problem (satellites, users, demands, etc.)
    problem = get_problem(config["date"], config["duration_factor"])

    # Solve the problem using all solvers
    results_all = {}
    for solver in solvers.keys():
        if solver not in config:
            continue
        results_all[solver] = solve_wrapper(config, problem, solver, n_jobs)
        logging.debug(f"\n{solver}", results_all[solver][0])

    return results_all


def get_problem(date, duration_factor) -> PntSchedulingProblem:
    """
    Create the PNT scheduling problem.

    Args:
        date (str): Start date
        duration_factor (float): Duration factor for demand generation

    Returns:
        PntSchedulingProblem: Problem
    """

    # *******************************************************************************
    # Satellite constellation
    # *******************************************************************************
    orbital_elements = pathfinder_data.orbital_elements.copy()
    users = pathfinder_data.users.copy()
    N_sat = 6
    N_planes = 3
    N_sat_plane = N_sat // N_planes

    # Classical orbital elements (a, e, i, W, w, M) [km, -, rad, rad, rad, rad]
    sma = 6142.4  # [km] Semi-major axis
    ecc = [0.001, 0.6, 0.6]  # [-] Eccentricity
    inc = np.deg2rad([0.001, 57.7, 57.7])  # [rad] Inclination
    raan = np.deg2rad([0, -90, 0])  # [rad] Right ascension of the ascending node
    aop = np.deg2rad([-45, -90, 90])  # [rad] Argument of periapsis
    ma = np.deg2rad([0, 180])  # [rad] Mean anomaly

    # Epoch (TAI)
    epoch_0 = pnt.SpiceInterface.string_to_tai(date)

    coe_OP = np.zeros((N_sat, 6))
    for i_pl in range(N_planes):
        for i_spl in range(N_sat_plane):
            coe_OP = np.array(
                [sma, ecc[i_pl], inc[i_pl], raan[i_pl], aop[i_pl], ma[i_spl]]
            )

    rv0_moon_sat_OP = pnt.classical_to_cartesian(coe_OP, pnt.MU_MOON)
    rv0_moon_sat_mi = pnt.CoordConverter.convert(
        epoch_0, rv0_moon_sat_OP, pnt.OP, pnt.MI
    )

    # Time
    sma = coe_OP[0, 0]  # [km] Semi-major axis
    period = 2 * np.pi * np.sqrt(np.power(sma, 3) / pnt.MU_MOON)  # [s] Orbital period
    Dt = 5 * pnt.SECS_PER_MINUTE  # [s] Simulation time step
    dt = 5 * pnt.SECS_PER_MINUTE  # [s] Propagation time step
    N_days = 2  # [-] Number of days
    tf = N_days * pnt.SECS_PER_DAY  # [s] Simulation final time
    N_t = int(tf / Dt)  # [-] Number of time steps
    tspan = np.linspace(0, tf, N_t)  # [s] Time since first epoch
    epochs = epoch_0 + tspan  # [s] Epochs (TAI)

    # Dynamics (three-body)
    dynamics = pnt.NBodyDynamics()
    dynamics.set_primary_body(pnt.Body.Moon())
    dynamics.add_body(pnt.Body.Earth())
    dynamics.set_time_step(dt)

    # Propagation
    # rv_{from}_{to}_{frame} [km, km/s] (x, y, z, vx, vy, vz)
    rv_moon_sat_mi = np.zeros((N_sat, N_t, 6))
    rv_moon_sat_pa = np.zeros((N_sat, N_t, 6))
    for i_sat in range(N_sat):
        rv_moon_sat_mi[i_sat] = dynamics.propagate(
            rv0_moon_sat_mi[i_sat], epoch_0, epochs
        )
        rv_moon_sat_pa[i_sat] = pnt.CoordConverter.convert(
            epochs, rv_moon_sat_mi[i_sat], pnt.MI, pnt.PA
        )
    rv_moon_earth_mi = pnt.SpiceInterface.get_body_pos_vel(epochs, pnt.MOON, pnt.EARTH)
    rv_moon_earth_pa = pnt.CoordConverter.convert(
        epochs, rv_moon_earth_mi, pnt.MI, pnt.PA
    )
    rv_moon_sun_mi = pnt.SpiceInterface.get_body_pos_vel(epochs, pnt.MOON, pnt.SUN)
    rv_moon_sun_pa = pnt.CoordConverter.convert(epochs, rv_moon_sun_mi, pnt.MI, pnt.PA)

    # Attitude
    r_sun = rv_moon_sun_mi[None, :, 0:3] - rv_moon_sat_mi[:, :, 0:3]
    r_earth = rv_moon_earth_mi[None, :, 0:3] - rv_moon_sat_mi[:, :, 0:3]

    e_sun = r_sun / np.linalg.norm(r_sun, axis=2)[:, :, None]
    e_earth = r_earth / np.linalg.norm(r_earth, axis=2)[:, :, None]

    r = rv_moon_sat_mi[:, :, 0:3]
    v = rv_moon_sat_mi[:, :, 3:6]
    r_norm = r / np.linalg.norm(r, axis=2)[:, :, None]
    v_norm = v / np.linalg.norm(v, axis=2)[:, :, None]
    e_nadir = -r_norm

    # Yaw-Steering (YS)
    ez_ys = -r_norm
    ey_ys = util.cross_norm(e_nadir, e_sun)
    ex_ys = util.cross_norm(ey_ys, ez_ys)

    # Sun angle
    sun_angle_cos = np.sum(e_sun * ex_ys, axis=2)
    sun_angle = np.arccos(sun_angle_cos)

    # *******************************************************************************
    # Users
    # *******************************************************************************

    def propagate_orbital_user(user: dict) -> np.array:
        """
        Propagate a user on orbit around the Moon.
        :param coe: [a, e, i, W, w, M] [km, -, deg, deg, deg, deg]
        :param frame: Frame
        :return: [x, y, z, vx, vy, vz] [km, km/s]
        """
        coe = user["orbital_elements"].copy()
        frame = user["frame"]
        coe[2:] = np.deg2rad(coe[2:])
        rv0 = pnt.classical_to_cartesian(coe, pnt.MU_MOON)
        rv0_mi = pnt.CoordConverter.convert(epoch_0, rv0, frame, pnt.MI)
        rv_mi = dynamics.propagate(rv0_mi, epoch_0, epochs)
        return rv_mi

    def propagate_surface_user(user: dict) -> np.array:
        """
        Propagate a user on the surface of the Moon.
        :param lat_lon_alt: [lat, lon, alt] [deg, deg, km]
        :return: [x, y, z, vx, vy, vz] [km, km/s]
        """
        lat_lon_alt = user["location"].copy()
        lat_lon_alt[:2] = np.deg2rad(lat_lon_alt[:2])
        rv_pa = np.zeros((N_t, 6))
        rv_pa[:, :3] = pnt.geographical_to_cartesian(lat_lon_alt, pnt.R_MOON)
        return rv_pa

    # Users
    N_users = len(users)
    user_type = np.array([user["type"] for user in users])
    rv_moon_user_mi = np.zeros((N_users, N_t, 6))
    rv_moon_user_pa = np.zeros((N_users, N_t, 6))
    az_el_rho = np.zeros((N_sat, N_users, N_t, 3))

    for i_usr, user in enumerate(users):
        if user["type"] == "orbital":
            rv_moon_user_mi[i_usr] = propagate_orbital_user(user)
            rv_moon_user_pa[i_usr] = pnt.CoordConverter.convert(
                epochs, rv_moon_user_mi[i_usr], pnt.MI, pnt.PA
            )
        elif user["type"] == "surface":
            rv_moon_user_pa[i_usr] = propagate_surface_user(user)
            rv_moon_user_mi[i_usr] = pnt.CoordConverter.convert(
                epochs, rv_moon_user_pa[i_usr], pnt.PA, pnt.MI
            )
        else:
            raise ValueError("Invalid user type")

        for j_sat in range(N_sat):
            az_el_rho[j_sat, i_usr] = pnt.cartesian_to_azimuth_elevation_range(
                rv_moon_user_mi[i_usr, :, :3], rv_moon_sat_mi[j_sat, :, :3]
            )

    # Elevation mask
    surface_elev_mask = np.deg2rad(15)  # [rad] Elevation mask
    orbital_elev_mask = np.deg2rad(0)  # [rad] Elevation mask
    max_elevation = np.max(az_el_rho[:, :, :, 1], axis=2)
    # min_elevation = np.maximum(surface_elev_mask, max_elevation - np.deg2rad(90))
    min_elevation = np.ones_like(max_elevation) * surface_elev_mask
    min_elevation[:, user_type == "orbital"] = orbital_elev_mask
    min_elevation[:, -1] = 0  # North Pole
    user_visible = np.greater_equal(az_el_rho[:, :, :, 1], min_elevation[:, :, None])

    # *******************************************************************************
    # Contacts
    # *******************************************************************************

    # Contact durations
    contact_durations = list[np.array]([] for _ in range(N_sat))
    contact_start_ends = list[np.array]([] for _ in range(N_sat))
    for i_sat in range(N_sat):
        for i, user in enumerate(users):
            starts, ends = util.get_start_end_indexes(user_visible[i_sat, i])
            contact_start_ends[i_sat].append(np.vstack((starts, ends)).T)
            contact_durations[i_sat].append((ends - starts) * Dt / pnt.SECS_PER_MINUTE)
    total_contact_durations = np.array(
        [[x.sum() for x in contact_durations[i_sat]] for i_sat in range(N_sat)]
    )
    contact_durations_pathfinder = np.array([user["contact"] for user in users])
    contact_durations_pathfinder = total_contact_durations[0] * duration_factor
    contact_durations_pathfinder = (
        np.round(contact_durations_pathfinder / 60 * 2) * 60 / 2
    )
    logging.debug(
        "Contact durations [minutes]",
        np.round(total_contact_durations / 60, 1),
        contact_durations_pathfinder / 60,
    )

    # Navigation signal
    fs = 2492.028e6  # Carrier frequency [Hz]
    fc = 5.115e6  # Spread code frequency [Hz]
    Tc = 1 / fc  # Spread code period [s]

    # Receiver parameters
    NF_lna = 1  # LNA noise figure [dB]
    T_sys = 113  # System noise temperature [K]
    B_dll = 0.5  # DLL bandwidth [Hz]
    d = 1.0  # Early-late spacing [chips]
    T_i = 0.02  # Coherent integration time [s]
    B_fe = 2 * fc  # Front-end bandwidth [Hz]

    # Constants
    c = 2.998e8  # Speed of light [m/s]
    k = -228.6  # Boltzmann constant [dBW/Hz/K]

    # Link budget
    rv_user_sat_mi = rv_moon_user_mi[None, :, :, :] - rv_moon_sat_mi[:, None, :, :]
    r = np.linalg.norm(rv_user_sat_mi[:, :, :, :3], axis=3)
    r[~user_visible] = np.nan
    # r = 8e6 # Distance [m]
    EIRP = 30  # Equivalent isotropic radiated power [dBW]
    L_fs = np.zeros_like(r)  # Free space loss [dB]
    for i_sat in range(N_sat):
        L_fs[i_sat] = pnt.decimal2dB(
            (4 * np.pi * r[i_sat] * fs / c) ** 2
        )  # Free space loss [dB]
    L_fs = 10 * np.log10((4 * np.pi * r * fs / c) ** 2)  # Free space loss [dB]
    Pr = EIRP - L_fs  # Received power [dBW]

    # Antenna gain
    T_atm = 290  # Atmospheric noise temperature [K]
    T_eq = 10 * np.log10(
        T_sys + T_atm * (10 ** (NF_lna / 10) - 1)
    )  # Equivalent noise temperature [dBK]

    mu = 0.6  # Antenna efficiency
    D_ant = 0.3  # Antenna diameter [m]
    # Maximum antenna gain [dBi]
    Gr_max = pnt.decimal2dB(mu * (np.pi * D_ant * fs / c) ** 2)
    Gr_max = 6

    theta = 0  # Beamwidth [deg]
    theta_3dB = 40  # 3 dB beamwidth [deg]
    Gr = pnt.decimal2dB(pnt.dB2decimal(Gr_max) - 12 * (theta / theta_3dB) ** 2)  # [dBi]

    # Gain-to-noise temperature [dB]
    GT = Gr - T_eq
    # GT = -20
    CN0 = Pr - k + GT  # Carrier-to-noise density ratio [dB-Hz]

    # Error contributions
    # Delay Lock Loop (DLL) error [m]
    term_1 = (c * Tc) ** 2 * (B_dll * (1 - 0.5 * B_dll * T_i)) / (2 * CN0)
    if d >= np.pi / (Tc * B_fe):
        sigma_p_dll_sq = term_1 * d
    elif np.pi / (Tc * B_fe) > d > 1 / (Tc * B_fe):
        term_2 = (1 / (Tc * B_fe)) + ((Tc * B_fe) / (np.pi - 1)) * (
            d - (1 / (Tc * B_fe))
        ) ** 2
        term_3 = 1 + (2 / (T_i * CN0 * (2 - d)))
        sigma_p_dll_sq = term_1 * term_2 * term_3
    else:  # d <= 1 / (T_c * B_fe)
        term_2 = 1 / (Tc * B_fe)
        term_3 = 1 + (1 / (T_i * CN0))
        sigma_p_dll_sq = term_1 * term_2 * term_3
    sigma_p_dll = np.sqrt(sigma_p_dll_sq)

    sigma_p_eph = 5  # Satellite ephemeris error [m]
    sigma_p_rel = 0.31  # Residual relay delay error [m]
    sigma_p_mul = 0.2  # Lunar multipath error [m]
    sigma_p_non_eph = np.sqrt(sigma_p_dll**2 + sigma_p_rel**2 + sigma_p_mul**2)
    sigma_p_tot = np.sqrt(4 * sigma_p_eph**2 + sigma_p_non_eph**2)

    # *******************************************************************************
    # Problem
    # *******************************************************************************

    # Requests
    requests = list[Request]()
    requests.append(
        Request(
            id=-1,
            usr_id=-1,
            rv=None,
            ts=0,
            te=tf / pnt.SECS_PER_HOUR,
            dur=tf / pnt.SECS_PER_HOUR,
            p=0,
        ),  # Dummy request
    )
    request_id = 0
    for j_day in range(N_days):
        for i, user in enumerate(users):
            requests.append(
                Request(
                    id=request_id,
                    usr_id=user["id"],
                    rv=rv_moon_user_mi[i],
                    ts=j_day * 24,
                    te=N_days * 24,
                    dur=contact_durations_pathfinder[i] / 2 / 60,
                )
            )
            request_id += 1
    N_req = len(requests)

    # Service windows
    service_windows = list[ServiceWindow]()
    window_id = 0
    for i_sat in range(N_sat):
        service_windows.append(
            ServiceWindow(
                id=-1,
                sat_id=i_sat,
                request_id=-1,
                ts=0,
                te=tf / pnt.SECS_PER_HOUR,
            )  # Dummy service window
        )
        for request in requests[1:]:
            for start, end in contact_start_ends[i_sat][request.usr_id]:

                # Convert indexes to hours
                start_h = start * Dt / pnt.SECS_PER_HOUR
                end_h = end * Dt / pnt.SECS_PER_HOUR

                if end_h <= request.ts or start_h >= request.te:
                    continue

                service_windows.append(
                    ServiceWindow(
                        id=window_id,
                        sat_id=i_sat,
                        request_id=request.id,
                        ts=ceil(start_h * 10) / 10,
                        te=floor(end_h * 10) / 10,
                    )
                )
                window_id += 1
    N_win = len(service_windows)

    # Transition times
    transition_times = np.ones((N_sat, N_win, N_win))
    transition_times *= 0.2  # hours
    for i_sat in range(N_sat):
        np.fill_diagonal(transition_times[i_sat, :, :], 0)
        transition_times[i_sat, 0, :] = 0
        transition_times[i_sat, :, 0] = 0

    # Energy and data generation functions
    S_panels = 3 * 0.5 * 0.8  # [m^2] Solar panels area
    I_earth = 1361  # [W/m^2] Solar irradiance
    efficiency = 0.2  # Solar panels efficiency
    P_solar_max = round(S_panels * I_earth * efficiency)  # [W] Max power generation
    P_solar_max = 50  # [W] Max power generation
    R_dte_max = -2e3  # [kbps] Direct-To-Earth max data rate

    payload_data_gen = -1.9 * R_dte_max  # [kbps] Payload data generation
    payload_energy_gen = -1.9 * P_solar_max  # [W] Payload power generation

    max_energy = np.abs(P_solar_max) * 12  # [Wh] Max energy capacity
    min_energy = max_energy * 0.2  # [Wh] Min energy capacity
    max_data = np.abs(R_dte_max) * 12  # [Mb] Max data capacity
    min_data = max_data * 0.2  # [Mb] Min data capacity

    e_func = partial(
        energy_gen_func,
        sun_angle_cos=sun_angle_cos,
        P_solar_max=P_solar_max,
        Dt=Dt,
    )
    d_func = partial(data_gen_func, R_dte_max=R_dte_max)

    # Problem
    problem = PntSchedulingProblem(
        time_step=Dt / pnt.SECS_PER_HOUR,
        requests=requests,
        service_windows=service_windows,
        transition_times=transition_times,
        CN0=CN0,
        max_energy=max_energy,
        min_energy=min_energy,
        max_data=max_data,
        min_data=min_data,
        payload_data_gen=payload_data_gen,
        payload_energy_gen=payload_energy_gen,
        energy_gen_func=e_func,
        data_gen_func=d_func,
    )

    return problem


def energy_gen_func(ts, te, sun_angle_cos, P_solar_max, Dt):
    """
    Satellite energy generation function.

    Args:
        ts (float): Start time [h]
        te (float): End time [h]
        sun_angle_cos (np.array): Sun angle cosine
        P_solar_max (float): Max power generation [W]
        Dt (float): Time step [s]

    Returns:
        float: Energy generation [Wh]
    """
    if not isinstance(ts, np.ndarray):
        i_s = round(ts * pnt.SECS_PER_HOUR / Dt)
        i_e = round(te * pnt.SECS_PER_HOUR / Dt)
        if i_e == i_s:
            if i_e < len(sun_angle_cos):
                i_e += 1
            else:
                i_s -= 1
        return P_solar_max * np.mean(sun_angle_cos[0, i_s:i_e]) * (te - ts)
    else:
        energy = np.zeros(len(ts))
        for i in range(len(ts)):
            i_s = round(ts[i] * pnt.SECS_PER_HOUR / Dt)
            i_e = round(te[i] * pnt.SECS_PER_HOUR / Dt)
            if i_e == i_s:
                if i_e < len(sun_angle_cos):
                    i_e += 1
                else:
                    i_s -= 1
            energy[i] = (
                P_solar_max * np.mean(sun_angle_cos[0, i_s:i_e]) * (te[i] - ts[i])
            )
        return energy


def data_gen_func(ts, te, R_dte_max):
    """
    Satellite data generation function (downlink).

    Args:
        ts (float): Start time [h]
        te (float): End time [h]
        R_dte_max (float): Max data rate [kbps]

    Returns:
        float: Data generation [Mb]
    """
    return R_dte_max * (te - ts)