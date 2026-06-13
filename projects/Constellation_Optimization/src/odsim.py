import numpy as np
import pylupnt as pnt
import os
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.optimize import curve_fit


class GroundStation:
    def __init__(self):
        """
        Initialize a ground station with latitude, longitude, and altitude.

        Parameters:
        lat (float): Latitude of the ground station in degrees.
        lon (float): Longitude of the ground station in degrees.
        alt (float): Altitude of the ground station in meters.
        """
        basepath = pnt.get_basepath()
        csv_file = os.path.join(basepath, "ground_station", "ground_stations.csv")
        df = pd.read_csv(csv_file)
        self.gs_ids = df["gs_id"].values
        self.latitudes = df["latitude_deg"].values
        self.longitudes = df["longitude_deg"].values
        self.altitudes = df["altitude_m"].values

    def get_latlonalt(self, gs_id):
        """
        Get the latitude, longitude, and altitude of the ground station by its ID.

        Parameters:
        gs_id (str): The ID of the ground station.

        Returns:
        tuple: (latitude in degrees, longitude in degrees, altitude in meters)
        """
        idx = np.where(self.gs_ids == gs_id)[0]
        if len(idx) == 0:
            raise ValueError(f"Ground station ID {gs_id} not found.")

        lat = self.latitudes[idx[0]]
        lon = self.longitudes[idx[0]]
        alt = self.altitudes[idx[0]]

        return np.array([lat, lon, alt])

    def get_ecef(self, gs_id):
        """
        Get the ECEF coordinates of the ground station by its ID.

        Parameters:
        gs_id (str): The ID of the ground station.

        Returns:
        tuple: ECEF coordinates (x, y, z) in meters.
        """
        idx = np.where(self.gs_ids == gs_id)[0]
        if len(idx) == 0:
            raise ValueError(f"Ground station ID {gs_id} not found.")

        lat_deg = self.latitudes[idx[0]]
        lon_deg = self.longitudes[idx[0]]
        alt_m = self.altitudes[idx[0]]

        RE = 6371.0  # Mean radius of the Earth in kilometers

        lat_rad = np.deg2rad(lat_deg)
        lon_rad = np.deg2rad(lon_deg)
        r_km = RE + alt_m / 1000.0  # Convert altitude to kilometers

        x = r_km * np.cos(lat_rad) * np.cos(lon_rad)
        y = r_km * np.cos(lat_rad) * np.sin(lon_rad)
        z = r_km * np.sin(lat_rad)

        pos = np.array([x, y, z])

        return pos


def propagate_sats_with_stm(rv0, et, dynamics=None):
    """
    Propagate satellite states with state transition matrix (STM).

    Parameters
    ----------
    rv0 : np.array (n_sat, 6)
        Initial satellite states in Moon-CI [r(3), v(3)].
    tspan : np.array
        Time steps for propagation.
    dynamics : pnt.Dynamics, optional
        Dynamics model for the simulation. If None, uses default NBodyDynamics.

    Returns
    -------
    x_sat : np.array (n_sat, T, 6)
        Propagated satellite states in Moon-CI.
    stm_sat : np.array (n_sat, T, 6, 6)
        State transition matrices for each satellite.
    """
    if dynamics is None:
        dynamics = pnt.NBodyDynamics()
        dynamics.set_integrator(pnt.IntegratorType.RKF45)
        dynamics.set_integrator_params(
            pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10)
        )
        dynamics.add_body(pnt.Body.Moon(20, 20))
        dynamics.add_body(pnt.Body.Earth())
        dynamics.add_body(pnt.Body.Sun())
        dynamics.set_time_step(60)  # 10 seconds time step
        dynamics.set_frame(pnt.MOON_CI)

    n_sat = rv0.shape[0]
    lent = et.shape[0]
    x_sat = np.zeros((n_sat, lent, 6))  # (N, T, 6) for [r(3), v(3)]
    stm_sat = np.zeros((n_sat, lent, 6, 6))  # (N, T, 6, 6) for state transition matrix

    # initialize the first state
    x_sat[:, 0, :] = rv0
    for i in range(n_sat):
        stm_sat[i, 0, :, :] = np.eye(6)  # identity matrix

    for i in range(1, lent):
        t0 = et[i - 1]
        tf = et[i]
        for j in range(n_sat):
            x_sat[j, i, :], stm_sat[j, i, :, :] = dynamics.propagate_stm(x_sat[j, i - 1, :], t0, tf)

    return x_sat, stm_sat


def compute_visibility(
    r1: np.ndarray, r2: np.ndarray, R_body: float, r_body: np.ndarray = None
) -> np.ndarray:
    if r_body is None:
        r_body = np.zeros(3)
    r = r2 - r1
    r_norm = np.linalg.norm(r, axis=1)
    r1body = r1 - r_body
    r1body_norm = np.linalg.norm(r1body, axis=1)
    dot = np.einsum("ij,ij->i", -r1body, r)
    theta1 = np.arccos(np.clip(dot / r_norm / r1body_norm, -1, 1))
    theta2 = np.arcsin(np.clip(R_body / r1body_norm, -1, 1))
    visibility = np.ones(len(r1), dtype=bool)
    visibility[(theta1 < theta2) & (r_norm > r1body_norm)] = False
    return visibility


def get_range_and_rate(
    et, gs_pos, sat_posvel_mci, get_H=False, add_noise=False, sigma_range=0.0, sigma_rangerate=0.0
):
    """
    Compute one-way range and range-rate between ground stations and satellites in ECEF.

    Parameters
    ----------
    et : np.ndarray
        Time steps in TAI.
    gs_pos : (N, 3) array
        Ground station positions [km] in ECEF.
    sat_posvel : (M, T, 6) array
        Satellite states [r (km), v (km/s)] in MCI,
    get_H : bool
        If True, also return Jacobian w.r.t. satellite state [r, v].
    add_noise : bool
        If True, add zero-mean Gaussian noise to the outputs.
    sigma_range : float
        Std dev for range noise [km].
    sigma_rangerate : float
        Std dev for range-rate noise [km/s].

    Returns
    -------
    elev: (N, M, T) array
        Visibility flag [0,1] for each GS/SAT pair.
    y : (2N, M, T) array
        First N rows: range [m]; next N rows: range-rate [km/s].
    H : (2N, M, T, 6) array, optional
        Jacobian of measurements w.r.t. satellite state [r(3), v(3)].
        Returned only if get_H is True.

    Notes
    -----
    - Uses one-way geometry: ρ = ||r_sat - r_gs||,  ρdot = u·v_sat (assuming v_gs = 0 in ECEF).
    - If you need two-way range, multiply the range block by 2 after computing.
    """
    # constants[
    N = gs_pos.shape[0]
    M = sat_posvel_mci.shape[0]
    T = et.shape[0]

    # first convert to ECEF
    Mrot = np.zeros((T, 6, 6))
    trans_vec = np.zeros((T, 6))
    sat_posvel_ecef = np.zeros((M, T, 6))  # (M, T, 6) for [r(3), v(3)]
    gs_pos_mci = np.zeros((N, T, 3))  # (N, 3)

    for i, t_tai in enumerate(et):
        Rot, rvec = pnt.get_frame_rotation_translation_rv(t_tai, pnt.MOON_CI, pnt.ECEF)
        sat_posvel_ecef[:, i, :] = (
            Rot @ sat_posvel_mci[:, i, :].transpose() + np.tile(rvec.reshape(6, 1), (1, M))
        ).transpose()  # (6, 6) x (6, M) + (6, M)
        gs_pos_mci[:, i, :] = pnt.convert_frame(t_tai, gs_pos, pnt.ECEF, pnt.MOON_CI)  # (N, 3)
        Mrot[i] = Rot
        trans_vec[i] = rvec

    sat_pos = sat_posvel_ecef[:, :, :3]  # (M, T, 3)
    sat_vel = sat_posvel_ecef[:, :, 3:]  # (M, T, 3)

    # compute the elevation
    # rho_vec: N x M x T x 3  (sat - gs)
    gs2e_vec = -gs_pos[:, np.newaxis, np.newaxis, :]  # N x 1 x 1 x 3
    gs2e_norm = np.tile(np.linalg.norm(gs2e_vec, axis=-1), (1, M, T))  # N x M x T
    rho_vec = sat_pos[np.newaxis, :, :, :] - gs_pos[:, np.newaxis, np.newaxis, :]  # N x M x T x 3
    rho = np.linalg.norm(rho_vec, axis=-1)
    dot_prod = np.sum(rho_vec * np.tile(gs2e_vec, (1, M, T, 1)), axis=-1)  # N x M x T
    elev = np.arccos(np.clip(dot_prod / rho / gs2e_norm, -1, 1)) - np.pi / 2  # N x M x T

    # moon occlusion
    vis_moon = np.ones((N, M, T), dtype=bool)
    for i in range(N):
        for j in range(M):
            vis_moon[i, j] = compute_visibility(
                gs_pos_mci[i, :, :3], sat_posvel_mci[j, :, :3], pnt.R_MOON
            )

    # Geometry
    # rho_vec: N x M x T x 3  (sat - gs)
    # N x M
    # Avoid divide-by-zero
    eps = 1e-12
    rho_safe = np.maximum(rho, eps)
    u = rho_vec / rho_safe[:, :, :, np.newaxis]  # N x M x T x 3

    # Range and range-rate (gs velocity assumed 0 in ECEF)
    rng = rho  # N x M x T
    v_rel = sat_vel[np.newaxis, :, :, :]  # N x M x T x 3
    rr = np.sum(u * v_rel, axis=-1)  # N x M x T

    # Optional noise
    if add_noise:
        if sigma_range > 0:
            rng = rng + np.random.normal(0.0, sigma_range, size=rng.shape)
        if sigma_rangerate > 0:
            rr = rr + np.random.normal(0.0, sigma_rangerate, size=rr.shape)

    # Stack outputs as (2N, M): first N rows = range, next N rows = range-rate
    y = np.vstack([rng, rr])

    if not get_H:
        return elev, y

    # Jacobians w.r.t. satellite state [r, v]  (per pair (i,j))
    # Range:  dρ/dr = u^T,  dρ/dv = 0
    H = np.zeros((2 * N, M, T, 6))
    H_range_pos = u  # N x M x T x 3
    H[:N, :, :, :3] = H_range_pos  # range w.r.t. position
    # H[:N, :, 3:] already zeros (range w.r.t. velocity)

    # Range-rate:  dρdot/dr = (1/ρ) * (I - u u^T) * v   (row vector)
    #              dρdot/dv = u^T
    # Compute dρdot/dr efficiently:
    # term = v - (u·v) u = v_rel - rr*u
    term = v_rel - rr[:, :, :, np.newaxis] * u  # N x M x T x 3
    H_rr_pos = term / rho_safe[:, :, :, np.newaxis]  # N x M x T x 3
    H_rr_vel = u  # N x M x T x 3

    H[N:, :, :, :3] = H_rr_pos
    H[N:, :, :, 3:] = H_rr_vel

    # Convert jacobian to with respect to satellite state [r(3), v(3)]
    for i in range(T):
        H[:, :, i, :] = (Mrot[i].T @ H[:, :, i, :].reshape(-1, 6).T).T.reshape(
            2 * N, M, 6
        )  # (6, 6) x (6, NxM) = (6, NxM)

    return elev, vis_moon, y, H


def Q_cv(sigma_a, dt):
    dt2 = dt * dt
    dt3 = dt2 * dt
    I3 = np.eye(3)
    Qpos = (dt3 / 3.0) * I3
    Qpv = (dt2 / 2.0) * I3
    Qvel = dt * I3
    top = np.hstack((Qpos, Qpv))
    bot = np.hstack((Qpv, Qvel))
    return (sigma_a**2) * np.vstack((top, bot))  # 6x6


def inertial_to_rtn_rotation(r, v):
    """
    Compute the 3x3 rotation matrix from inertial (ECI) to RTN frame.

    Parameters
    ----------
    state_inertial : array_like, shape (N, 6)
        Inertial state vector [rx, ry, rz, vx, vy, vz] in meters and m/s.

    Returns
    -------
    R_eci2rtn : ndarray, shape (N,3,3)
        Rotation matrix that transforms a vector in ECI to RTN frame.
        RTN axes definition:
            R: radial (along position vector)
            T: along-track (in orbital plane, orthogonal to R, in direction of velocity)
            N: cross-track (completes right-handed system)
    """
    # Unit radial vector
    r_hat = r / np.linalg.norm(r, axis=-1, keepdims=True)  # (N, 3)

    # Angular momentum vector (h = r × v)
    h = np.cross(r, v)  # (N, 3)
    h_norm = np.linalg.norm(h, axis=-1, keepdims=True)  # (N, 1)
    n_hat = h / h_norm  # Normalized angular momentum vector (N, 3)

    # Along-track unit vector
    t_hat = np.cross(n_hat, r_hat)  # (N, 3)

    # Build rotation matrix (columns are RTN basis in ECI)
    R_eci2rtn = np.stack((r_hat, t_hat, n_hat), axis=-1)  # (N, 3, 3)

    return R_eci2rtn


class ODSim:
    def __init__(self, od_config):
        """
        Initialize the ODSim with the given configuration.

        Parameters:
        od_config (dict): Configuration dictionary for the ODSim.
        """
        self.od_config = od_config
        self.et0 = od_config["et0"]
        self.meas = od_config["meas"]

        # setup ground stations
        if "gs" in self.meas:
            self.gs = GroundStation()
            self.n_gs = len(self.od_config["gs_ids"])
            self.pos_gs = np.zeros((self.n_gs, 3))
            self.lla_gs = np.zeros((self.n_gs, 3))
            for i, gs_id in enumerate(od_config["gs_ids"]):
                self.lla_gs[i] = self.gs.get_latlonalt(gs_id)
                self.pos_gs[i] = self.gs.get_ecef(gs_id)

    def run_sim(self, tspan, rv0=None, dynamics=None, res=None):
        """_summary_

        Parameters
        ----------
        tspan : np.array
            time steps
        rv0 : np.array (n_sat, 6)
            satellite states in Moon-CI
        dyn : pnt.Dynamics
            dynamics model for the simulation
        """
        et = self.od_config["et0"] + tspan  # convert to TAI
        dt = tspan[1] - tspan[0]
        len_prop = int(
            self.od_config["predict_time"] / dt
        )  # number of steps to propagate after measurement update

        # first, propagate the satellite states ------------------------------------
        if res is None:
            if rv0 is None:
                raise ValueError("rv0 must be provided if res is None")
            x_sat, stm_sat = propagate_sats_with_stm(rv0, et, dynamics=dynamics)

            # generate all the Jacobians throught the trajectory -----------------------------------
            elev, vis_moon, y, H = get_range_and_rate(
                et,
                self.pos_gs,
                x_sat,
                get_H=True,
                add_noise=self.od_config["add_noise"],
                sigma_range=self.od_config["sigma_range"],
                sigma_rangerate=self.od_config["sigma_rangerate"],
            )
        else:
            x_sat = res["x_sat"]
            coe = res["coe"]
            stm_sat = res["stm_sat"]
            elev = res["elev"]
            vis_moon = res["vis_moon"]
            y = res["y"]
            H = res["H"]
            rv0 = x_sat[:, 0, :]  # initial state (n_sat, 6)

        # constants
        n_sat = rv0.shape[0]
        n_gs = self.n_gs
        lent = tspan.shape[0]

        # set initial covariance ----------------------------------------------------------
        P0 = np.diag(
            np.hstack(
                [
                    self.od_config["init_pos_std"] ** 2 * np.ones(3),
                    self.od_config["init_vel_std"] ** 2 * np.ones(3),
                ]
            )
        )

        Psats = np.zeros((n_sat, 6, 6))  # (N, 6, 6) for each satellite
        for i in range(n_sat):
            Psats[i] = P0

        dt = tspan[1] - tspan[0]  # time step in seconds

        # SRP perturbation model
        SOLAR_FLUX_AU = 1367.0  # [W/m^2] Solar constant at 1 AU
        C_MS = 299792458.0  # [m/s] Speed of light
        P_SUN = SOLAR_FLUX_AU / C_MS  # [N/m^2] Solar power per unit area at 1 AU
        CR = 1.8  # [-] Coefficient of reflectivity (assumed constant)
        area = 1.0  # [m^2] Area of the satellite
        mass = 850.0  # [kg] Mass of the satellite
        coeff = CR * area / mass  # [N/kg = m/s^2] SRP acceleration coefficient
        srp = coeff * P_SUN  # [N/kg = m/s^2] SRP acceleration coefficient

        # Q = Q_cv(srp*1e-3, dt)  # process noise
        Q = Q_cv(self.od_config["std_acc"] * 1e-3, dt)  # process noise

        # MCI -> RTN rotation matrix
        R_eci2rtn = inertial_to_rtn_rotation(x_sat[:, :, :3], x_sat[:, :, 3:6])  # (M, T, 3)

        # then run covariance analysis ----------------------------------------------------------
        sigma_sat_mci = np.zeros((n_sat, lent, 6))  # (N, T, 6) for each satellite
        sigma_sat_rtn = np.zeros((n_sat, lent, 6))

        for i, t_tai in enumerate(tspan):
            elev_t = elev[:, :, i]
            vis_t = (elev_t >= self.od_config["elev_mask_deg"] * np.pi / 180.0).astype(
                int
            )  # elevation mask
            vis_moon_t = vis_moon[:, :, i]
            vis_t = vis_t * vis_moon_t  # (N, M) final visibility mask
            vis_t = np.tile(vis_t, (2, 1))  # for range and range-rate

            # extract the measurement at this time step
            y_t = y[:, :, i]  # (2N, M)
            H_t = H[:, :, i, :]  # (2N, M, 6)

            for j in range(n_sat):
                Psat = Psats[j]

                # Dynamics update ----------------------------------------------------------
                # propagate the covariance using the STM
                stm = stm_sat[j, i, :, :]  # (6, 6)
                Psat = stm @ Psat @ stm.T + Q

                # Measurement update ----------------------------------------------------------
                y_vis = y_t[:, j][vis_t[:, j] == 1]
                H_vis = H_t[:, j, :][vis_t[:, j] == 1, :]
                n_vis = int(y_vis.shape[0] / 2)

                if n_vis > 0:
                    # compute the measurement covariance
                    R = np.diag(
                        np.hstack(
                            [
                                self.od_config["sigma_range"] ** 2 * np.ones(n_vis),
                                self.od_config["sigma_rangerate"] ** 2 * np.ones(n_vis),
                            ]
                        )
                    )

                    # compute the Kalman gain
                    S = H_vis @ Psat @ H_vis.T + R
                    K = Psat @ H_vis.T @ np.linalg.pinv(S)

                    # Update covariance in Jacobian form
                    G = np.eye(6) - K @ H_vis
                    Psat = G @ Psat @ G.T + K @ R @ K.T

                Psats[j] = Psat

                # propagate the covariance to predict future steps
                Psat_prop = Psat.copy()
                for k in range(i + 1, i + len_prop):
                    if k >= lent:
                        break
                    stm = stm_sat[j, k, :, :]
                    Psat_prop = stm @ Psat_prop @ stm.T + Q

                # convert to RTN frame
                i2rtn = R_eci2rtn[j, i, :, :]
                i2rtn = np.block([[i2rtn, np.zeros((3, 3))], [np.zeros((3, 3)), i2rtn]])
                Psat_rtn = i2rtn @ Psat_prop @ i2rtn.T

                # store (store the worst-case covariance before measurement update)
                sigma_sat_mci[j, i, :] = np.sqrt(np.diag(Psat_prop))
                sigma_sat_rtn[j, i, :] = np.sqrt(np.diag(Psat_rtn))

        return sigma_sat_mci, sigma_sat_rtn


def rms(x, axis=None):
    """Calculate the root mean square of a vector."""
    return np.sqrt(np.mean(x**2, axis=axis))


def plot_od_result(
    tspan,
    elev,
    vis_moon,
    sigma_sat_mci,
    sigma_sat_rtn,
    elev_mask_deg=10,
    stats_ratio=2 / 3,
    use_logy=False,
    figname=None,
):

    n_sat = sigma_sat_mci.shape[0]
    lent = tspan.shape[0]

    fig, ax = plt.subplots(n_sat, 2, figsize=(10, 3 * n_sat), sharex=False)

    if n_sat == 1:
        ax = ax[np.newaxis, :]

    print("Sat        R [m]          T [m]          N [m]       Norm [m]")
    print("---------------------------------------------------------------")
    for i in range(n_sat):
        vis_t = elev[:, i, :] >= elev_mask_deg * np.pi / 180.0
        vis_t_moon = vis_moon[:, i, :]
        vis_t = vis_t * vis_t_moon
        vis_num = np.sum(vis_t, axis=0)
        ax[i, 0].plot(tspan / 3600, 1000 * sigma_sat_rtn[i, :, 0], label="R")
        ax[i, 0].plot(tspan / 3600, 1000 * sigma_sat_rtn[i, :, 1], label="T")
        ax[i, 0].plot(tspan / 3600, 1000 * sigma_sat_rtn[i, :, 2], label="N")
        ax[i, 0].plot(
            tspan / 3600,
            1000 * np.linalg.norm(sigma_sat_mci[i, :, :3], axis=1),
            label="Norm",
            linestyle="--",
            color="black",
        )
        # gray out the times when no measurement is available
        ax[i, 0].fill_between(
            tspan / 3600,
            0,
            100,
            where=vis_num == 0,
            color="gray",
            alpha=0.3,
            label="No Measurement",
        )
        ax[i, 0].set_ylabel(f"RTN Position Error (m)", fontsize=12)
        ax[i, 0].set_xlabel("Time (hours)", fontsize=12)
        ax[i, 0].grid()
        ax[i, 0].legend()

        ax[i, 1].plot(tspan / 3600, 1e6 * sigma_sat_rtn[i, :, 3], label="V_R")
        ax[i, 1].plot(tspan / 3600, 1e6 * sigma_sat_rtn[i, :, 4], label="V_T")
        ax[i, 1].plot(tspan / 3600, 1e6 * sigma_sat_rtn[i, :, 5], label="V_N")
        ax[i, 1].plot(
            tspan / 3600,
            1e6 * np.linalg.norm(sigma_sat_mci[i, :, 3:], axis=1),
            label="Norm",
            linestyle="--",
            color="black",
        )
        ax[i, 1].fill_between(
            tspan / 3600,
            0,
            100,
            where=vis_num == 0,
            color="gray",
            alpha=0.3,
            label="No Measurement",
        )
        ax[i, 1].set_ylabel(f"RTN Velocity Error (mm/s)", fontsize=12)
        ax[i, 1].set_xlabel("Time (hours)", fontsize=12)
        ax[i, 1].grid()
        ax[i, 1].legend()

        if use_logy:
            ax[i, 0].set_yscale("log")
            ax[i, 1].set_yscale("log")
        else:
            ax[i, 0].set_ylim((0, 50))
            ax[i, 1].set_ylim((0, 5))

        # statistics for last 1/3 of the data
        start_idx = int(lent * stats_ratio)
        pos_rtn_rms = rms(sigma_sat_rtn[i, start_idx:, :], axis=0) * 1000
        pos_rtn_95 = np.percentile(sigma_sat_rtn[i, start_idx:, :], 95, axis=0) * 1000
        pos_norm_rms = rms(np.linalg.norm(sigma_sat_mci[i, start_idx:, :3], axis=1)) * 1000
        pos_norm_95 = (
            np.percentile(np.linalg.norm(sigma_sat_mci[i, start_idx:, :3], axis=1), 95) * 1000
        )
        print(
            f"SAT {i+1}  | "
            f"{pos_rtn_rms[0]:.2f} ({pos_rtn_95[0]:.2f})   "
            f"{pos_rtn_rms[1]:.2f} ({pos_rtn_95[1]:.2f})  "
            f"{pos_rtn_rms[2]:.2f} ({pos_rtn_95[2]:.2f})   "
            f"{pos_norm_rms:.2f} ({pos_norm_95:.2f}) "
        )

    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname)

    plt.show()


def fit_od_error(tspan, T_orbit, sigma_sat_rtn, stats_ratio=2 / 3, debug=False, figname=None):
    """
    Fit the OD error over one orbit to a decaying exponential function.

    Parameters
    ----------
    tspan : np.array
        Time steps.
    T_orbit : float
        Orbital period in seconds.
    sigma_sat_rtn : np.array (n_sat, T, 6)
        RTN position/velocity errors for each satellite.
    stats_ratio : float
        Fraction of the orbit to consider for fitting.

    Returns
    -------
    fit_params : list of dict
        List of fitted parameters for each satellite.
    """
    n_params = 10  # number of parameters to fit

    def curve_func(x, a0, a1, c1, s1, c2, s2, c3, s3, c4, s4):
        y1 = a0 + a1 * x
        y2 = (
            s1 * np.sin(x)
            + c1 * np.cos(x)
            + s2 * np.sin(2 * x)
            + c2 * np.cos(2 * x)
            + s3 * np.sin(3 * x)
            + c3 * np.cos(3 * x)
            + c4 * np.cos(4 * x)
            + s4 * np.sin(4 * x)
        )
        return y1 + y2

    n_sat = sigma_sat_rtn.shape[0]
    lent = sigma_sat_rtn.shape[1]
    fit_params = np.zeros((n_sat, n_params))

    if debug:
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), sharex=True)

    start_idx = int(lent * stats_ratio)
    for i in range(n_sat):
        t_fit = (
            (tspan[start_idx:] - tspan[start_idx]) / T_orbit * 2 * np.pi
        )  # normalize to [0, 2pi]
        y_fit = sigma_sat_rtn[i, start_idx:, 0]

        # Fit each RTN component separately
        try:
            popt, _ = curve_fit(curve_func, t_fit, y_fit)
            params = popt
        except RuntimeError:
            params = None
            print(f"Warning: Fit failed for satellite {i+1}")

        if debug:
            ax.plot(t_fit, y_fit * 1000, "o-", label=f"Sat {i+1} R Data")
            if params is not None:
                y_model = curve_func(t_fit, *params)
                ax.plot(t_fit, 1000 * y_model, linestyle="--", label=f"Sat {i+1} R Fit")
            ax.set_ylabel("R Position Error (m)", fontsize=12)
            ax.set_xlabel("Time (hours)", fontsize=12)
            ax.grid()
            ax.legend()

        fit_params[i, :] = params if params is not None else np.nan

    if debug:
        plt.tight_layout()
        if figname is not None:
            plt.savefig(figname)
        plt.show()

    return fit_params


# --- Worker: computes metrics for one (sma, inc) cell ---
def compute_one_cell(i, j, index, sim_num, sma, inc, omegas, et0, n_orbit=4):
    """
    Return: (i, j, metrics) where metrics has shape (n_omega, 8)
    dyn_cfg / odsim_cfg should be PICKLABLE configs or factory args,
    not heavy unpicklable objects. Create the objects inside the worker.
    """
    # (Re)create dyn / odsim here from configs (recommended for pickling)
    odconfig = {
        "et0": et0,
        "meas": ["gs"],
        "gs_ids": ["LEGS1_X", "LEGS2_X", "DSS35_X"],
        "init_pos_std": 100e-3,  # km
        "init_vel_std": 1e-3,  # km/s
        "sigma_range": 1e-3,  # m
        "sigma_rangerate": 0.1e-6,  # km/s
        "elev_mask_deg": 10,  # degrees
        "std_acc": 1e-7,  # km/s^2
        "predict_time": 2 * 3600,  # seconds
    }

    odsim = ODSim(odconfig)
    print(
        "Running simulation {}/{} for sma={:.1f} km, inc={:.1f} deg".format(
            index + 1, sim_num, sma, np.rad2deg(inc)
        )
    )

    # dynamics setup
    dyn = pnt.NBodyDynamics()
    dyn.set_integrator(pnt.IntegratorType.RKF45)
    dyn.set_integrator_params(pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10))
    dyn.add_body(pnt.Body.Moon(20, 20))
    dyn.add_body(pnt.Body.Earth())
    dyn.add_body(pnt.Body.Sun())
    dyn.set_time_step(60)  # 10 seconds time step
    dyn.set_frame(pnt.MOON_CI)
    dyn.set_autodiff(True)

    n_sat = len(omegas)
    results_ij = np.ones((n_sat, 8), dtype=float) * np.nan

    # Skip infeasible eccentricities
    ecc = np.sqrt(max(0.0, 1 - 5 / 3 * np.cos(inc) ** 2))
    peri_h = sma * (1 - ecc) - pnt.R_MOON
    if peri_h < 100:
        return i, j, results_ij

    w = np.deg2rad(-90.0)
    M = np.deg2rad(0.0)

    T_orbit = 2 * np.pi * np.sqrt(sma**3 / pnt.GM_MOON)
    dt_od = 60.0
    n_od = int(T_orbit / dt_od)
    tspan = np.linspace(0.0, n_orbit * T_orbit, n_od + 1)
    et = et0 + tspan

    # COEs for all sats at this cell
    coes = np.zeros((n_sat, 6))
    coes[:, 0] = sma
    coes[:, 1] = ecc
    coes[:, 2] = inc
    coes[:, 3] = omegas
    coes[:, 4] = w
    coes[:, 5] = M

    # Initial states in MOON_CI
    rv0_mci = np.zeros((n_sat, 6))
    for si in range(n_sat):
        rv0_op = pnt.classical_to_cart(coes[si, :], pnt.GM_MOON)
        rv0_mci[si] = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)

    # Propagate all sats together
    x_sat, stm_sat = propagate_sats_with_stm(rv0_mci, et, dynamics=dyn)

    # Generate measurements
    elev, vis_moon, y, H = get_range_and_rate(
        et,
        odsim.pos_gs,
        x_sat,
        get_H=True,
        add_noise=True,
        sigma_range=1e-3,
        sigma_rangerate=1e-7,
    )

    # Run OD sim (returns 1-sigma cov-derived stds per axis)
    res = {
        "x_sat": x_sat,
        "stm_sat": stm_sat,
        "elev": elev,
        "vis_moon": vis_moon,
        "y": y,
        "H": H,
        "coe": coes,
    }
    sigma_sat_mci, sigma_sat_rtn = odsim.run_sim(tspan, rv0_mci, dynamics=dyn, res=res)

    # Compute metrics only over the latter fraction of the arc
    stats_ratio = (n_orbit - 1) / n_orbit
    start_idx = int(len(tspan) * stats_ratio)

    for k in range(n_sat):
        pos_rtn = sigma_sat_rtn[k, start_idx:, :3]  # (T, 3) [km]
        pos_mci = sigma_sat_mci[k, start_idx:, :3]  # (T, 3) [km]

        pos_rtn_rms = np.sqrt(np.mean(pos_rtn**2, axis=0)) * 1000.0  # m
        pos_norm_rms = np.sqrt(np.mean(np.linalg.norm(pos_mci, axis=1) ** 2)) * 1000.0  # m
        pos_rtn_95 = np.percentile(pos_rtn, 95, axis=0) * 1000.0  # m
        pos_norm_95 = np.percentile(np.linalg.norm(pos_mci, axis=1), 95) * 1000.0  # m

        metrics = np.empty(8, dtype=float)
        metrics[0:3] = pos_rtn_95
        metrics[3] = pos_norm_95
        metrics[4:7] = pos_rtn_rms
        metrics[7] = pos_norm_rms

        results_ij[k, :] = metrics

    return i, j, results_ij


def gridsearch_od_parallel(et0, smas, incs, omegas, n_orbit=4, max_workers=None):
    """
    Parallel version. odsim_cfg and dyn_cfg should be either:
      - picklable lightweight objects, or
      - zero-arg factory callables that build the heavy objects in the worker.
    """
    n_sma = len(smas)
    n_incs = len(incs)
    n_omega = len(omegas)

    sim_num = n_sma * n_incs

    results = np.ones((n_sma, n_incs, n_omega, 8), dtype=float) * np.nan

    tasks = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        for i, sma in enumerate(smas):
            for j, inc in enumerate(incs):
                tasks.append(
                    ex.submit(
                        compute_one_cell,
                        i,
                        j,
                        i * n_incs + j,
                        sim_num,
                        float(sma),
                        float(inc),
                        np.asarray(omegas, dtype=float),
                        float(et0),
                        int(n_orbit),
                    )
                )

        for fut in as_completed(tasks):
            i, j, metrics = fut.result()
            results[i, j, :, :] = metrics

    return results


def fit_params_one_cell(i, j, index, sim_num, sma, ecc, omegas, et0, n_orbit=4):
    """
    Return: (i, j, fit_params) where fit_params has shape (n_omega, 10)
    dyn_cfg / odsim_cfg should be PICKLABLE configs or factory args,
    not heavy unpicklable objects. Create the objects inside the worker.
    """
    # (Re)create dyn / odsim here from configs (recommended for pickling)
    odconfig = {
        "et0": et0,
        "meas": ["gs"],
        "gs_ids": ["LEGS1_X", "LEGS2_X", "DSS35_X"],
        "init_pos_std": 100e-3,  # km
        "init_vel_std": 1e-3,  # km/s
        "sigma_range": 1e-3,  # m
        "sigma_rangerate": 0.1e-6,  # km/s
        "elev_mask_deg": 10,  # degrees
        "std_acc": 1e-7,  # km/s^2
        "predict_time": 2 * 3600,  # seconds
    }

    odsim = ODSim(odconfig)
    print(
        "Running simulation {}/{} for sma={:.1f} km, ecc={:.2f}".format(
            index + 1, sim_num, sma, ecc
        )
    )

    # dynamics setup
    dyn = pnt.NBodyDynamics()
    dyn.set_integrator(pnt.IntegratorType.RKF45)
    dyn.set_integrator_params(pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10))
    dyn.add_body(pnt.Body.Moon(20, 20))
    dyn.add_body(pnt.Body.Earth())
    dyn.add_body(pnt.Body.Sun())
    dyn.set_time_step(60)  # 10 seconds time step
    dyn.set_frame(pnt.MOON_CI)
    dyn.set_autodiff(True)

    n_sat = len(omegas)
    results_ij = np.ones((n_sat, 10), dtype=float) * np.nan

    # Skip infeasible eccentricities
    peri_h = sma * (1 - ecc) - pnt.R_MOON
    if peri_h < 100:
        return i, j, results_ij  # return NaNs

    inc = np.arccos(np.sqrt(3 / 5 * (1 - ecc**2)))  # rad

    w = np.deg2rad(90.0)
    M = np.deg2rad(0.0)
    T_orbit = 2 * np.pi * np.sqrt(sma**3 / pnt.GM_MOON)
    dt_od = 60.0
    n_od = int(T_orbit / dt_od)
    tspan = np.linspace(0.0, n_orbit * T_orbit, n_od + 1)
    et = et0 + tspan
    # COEs for all sats at this cell
    coes = np.zeros((n_sat, 6))
    coes[:, 0] = sma
    coes[:, 1] = ecc
    coes[:, 2] = inc
    coes[:, 3] = omegas
    coes[:, 4] = w
    coes[:, 5] = M
    # Initial states in MOON_CI
    rv0_mci = np.zeros((n_sat, 6))
    for si in range(n_sat):
        rv0_op = pnt.classical_to_cart(coes[si, :], pnt.GM_MOON)
        rv0_mci[si] = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)
    # Propagate all sats together
    x_sat, stm_sat = propagate_sats_with_stm(rv0_mci, et, dynamics=dyn)
    # Generate measurements
    elev, vis_moon, y, H = get_range_and_rate(
        et,
        odsim.pos_gs,
        x_sat,
        get_H=True,
        add_noise=False,
        sigma_range=1e-3,
        sigma_rangerate=1e-7,
    )
    # Run OD sim (returns 1-sigma cov-derived stds per axis)
    res = {
        "x_sat": x_sat,
        "stm_sat": stm_sat,
        "elev": elev,
        "vis_moon": vis_moon,
        "y": y,
        "H": H,
        "coe": coes,
    }
    sigma_sat_mci, sigma_sat_rtn = odsim.run_sim(tspan, rv0_mci, dynamics=dyn, res=res)
    # Fit params only over the latter fraction of the arc
    stats_ratio = (n_orbit - 1) / n_orbit
    fit_params = fit_od_error(tspan, T_orbit, sigma_sat_rtn, stats_ratio=stats_ratio, debug=False)

    return i, j, fit_params


def fit_oderr_params_parallel(et0, smas, eccs, omegas, n_orbit=4, max_workers=None):
    """
    Parallel version. odsim_cfg and dyn_cfg should be either:
      - picklable lightweight objects, or
      - zero-arg factory callables that build the heavy objects in the worker.
    """
    n_sma = len(smas)
    n_eccs = len(eccs)
    n_omega = len(omegas)
    n_params = 10

    sim_num = n_sma * n_eccs

    fit_params_all = np.ones((n_sma, n_eccs, n_omega, n_params), dtype=float) * np.nan

    tasks = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        for i, sma in enumerate(smas):
            for j, ecc in enumerate(eccs):
                tasks.append(
                    ex.submit(
                        fit_params_one_cell,
                        i,
                        j,
                        i * n_eccs + j,
                        sim_num,
                        float(sma),
                        float(ecc),
                        np.asarray(omegas, dtype=float),
                        float(et0),
                        int(n_orbit),
                    )
                )

        for fut in as_completed(tasks):
            i, j, params = fut.result()
            fit_params_all[i, j, :, :] = params

    return fit_params_all


def param_to_oderr(x, param):
    """
    Compute the orbit determination error based on the fitting parameters.

    Parameters:
    t : np.ndarray
        Time since epoch in seconds.
    param : list or np.ndarray
        Fitting parameters [a0, a1, a2, b0, b1, b2].

    Returns:
    od_err : np.ndarray
        Orbit determination error in meters.
    """
    a0, a1, c1, s1, c2, s2, c3, s3, c4, s4 = param
    y1 = a0 + a1 * x
    y2 = (
        s1 * np.sin(x)
        + c1 * np.cos(x)
        + s2 * np.sin(2 * x)
        + c2 * np.cos(2 * x)
        + s3 * np.sin(3 * x)
        + c3 * np.cos(3 * x)
        + c4 * np.cos(4 * x)
        + s4 * np.sin(4 * x)
    )
    y = y1 + y2
    return y
