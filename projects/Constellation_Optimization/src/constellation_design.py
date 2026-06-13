import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pymoo.core.problem import ElementwiseProblem
from tqdm import tqdm
import time
import itertools
from src.failure_model import NavSatFailureModel
from src.communication import LunanetSatAntenna
from src.opt_objectives import compute_coverage_dop
from pymoo.core.variable import Real, Integer, Binary
import os

# constants
N_WALKER_PARAMS = 7  # number of parameters per walker constellation
N_WALKER_CONSTRAINTS = 2  # number of constraints per walker constellation
N_GLOBAL_CONSTRAINTS = 2  # number of global constraints
OBJ_INF = 10.0
N_SAT_MAX = 50
N_PHASE_VEC = 60  # phase vector
P_TX_W = 15.0  # [dBW] transmission power of each satellite
OPT_BASEPATH = "/Users/keidaiiiyama/Documents/sw_navlab/LuPNT-private/python/examples/ex_constellation_optimization/"
N_SMA_INT = 7  # number of discrete sma integer values (1/1, 2/3, 1/2, 2/5, 1/3, 2/7, 1/4) of 48 hr


class ConstellationOptimization(ElementwiseProblem):

    def __init__(self, config, verbose, **kwargs):
        self.config = config
        self.n_walker = config["n_walker"]
        self.n_obj = config["n_obj"]
        self.n_var = config["n_var"]
        self.verbose_ = verbose
        self.evalcount = 0

        self._set_vars()

        super().__init__(
            vars=self.vars,
            n_obj=self.config["n_obj"],
            n_ieq_constr=self.config["n_ieq_constr"],
            **kwargs,
        )

    def _evaluate(self, X, out, *args, **kwargs):
        """
        Evaluate the objective and constraints for a given design variable x

        Parameters:
        x: design vector [8 * n_walker]

        Returns:
        enu_mats : np.ndarray [n_users, 3, 3]
            ENU matrices for each user position.
        """
        x = self._get_xvector(X)

        out["G"], params_c = compute_constraints(x, self.config)

        if np.any(out["G"] > 0):
            # print("Constraints violated!")
            out["F"] = np.ones(self.n_obj) * OBJ_INF
        else:
            # print("Computing objectives...")
            out["F"], labels, params_o = compute_objectives(x, self.config)

        # print_lunanet_x(x, self.config)
        if self.verbose_:
            self.evalcount += 1
            print("eval:", self.evalcount, "  objectives: ", out["F"])

    def _set_vars(self):
        keys = []
        vars = {}
        xl = self.config["xl"]
        xu = self.config["xu"]
        float_sma = self.config.get("float_sma", True)  # optimize semi-major axis
        k = 0
        masks = []

        for i in range(self.n_walker):
            # semi-major axis (real)
            if float_sma:
                vars["a_{}".format(i)] = Real(bounds=(xl[k], xu[k]))
                keys.append("a_{}".format(i))
                masks.append("real")
            else:
                # instead optimize integer multiple of orbital period
                vars["a_{}".format(i)] = Integer(bounds=(xl[k], xu[k]))
                keys.append("a_{}".format(i))
                masks.append("int")
            k += 1
            # 1 eccentricity (real)
            vars["e_{}".format(i)] = Real(bounds=(xl[k], xu[k]))
            keys.append("e_{}".format(i))
            masks.append("real")
            k += 1
            # 2 sign of argument of periapsis (binary)
            vars["i_{}".format(i)] = Binary()
            keys.append("i_{}".format(i))
            masks.append("bool")
            k += 1
            # 3 number of planes (integer)
            vars["n_{}".format(i)] = Integer(bounds=(xl[k], xu[k]))
            keys.append("n_{}".format(i))
            masks.append("int")
            k += 1
            # 4 number of satellites per plane (integer)
            vars["m_{}".format(i)] = Integer(bounds=(xl[k], xu[k]))
            keys.append("m_{}".format(i))
            masks.append("int")
            k += 1
            # 5 relative phase (integer)
            vars["phi_{}".format(i)] = Integer(bounds=(xl[k], xu[k]))
            keys.append("phi_{}".format(i))
            masks.append("int")
            k += 1
            # 6 Omega 0 (real)
            vars["Omega_0_{}".format(i)] = Real(bounds=(xl[k], xu[k]))
            keys.append("Omega_0_{}".format(i))
            masks.append("real")
            k += 1

        # phase allocation variable (integer)
        for i in range(N_PHASE_VEC):
            vars["phase_{}".format(i)] = Integer(bounds=(-1, self.config["n_phase"] - 1))
            keys.append("phase_{}".format(i))
            masks.append("int")

        self.vars = vars
        self.keys = keys
        self.masks = masks

    def _get_xvector(self, X):
        x = np.array([X[key] for key in self.keys])
        # print_lunanet_x(x, self.config)
        return x

    def _x_from_vector(self, x_vector):
        X = {}
        for i, key in enumerate(self.keys):
            if self.masks[i] == "int":
                X[key] = int(x_vector[i])
            elif self.masks[i] == "bool":
                X[key] = bool(x_vector[i])
            else:
                X[key] = x_vector[i]
        return X


def setup_problem_config(
    objs,
    n_walker=3,
    n_phase=1,
    sat_range_phases=[[4, 10]],
    float_sma=True,
    n_users=100,
    elev_mask_deg=5.0,
    cn0_thresh_user=30.0,
    dop_type_phase=None,
    use_wdop=False,
    dop_target=None,
    compute_link_budget=True,
    lat_masks=[[-90, -75]],
    launch_years=[0],
    eval_years=[1],
    fail_model=NavSatFailureModel(),
    max_fail_sat=1,
    epoch=None,
    tspan=None,
    dt=None,
):
    """
    Set up the problem configuration for constellation optimization.

    Parameters
    --------

    """

    objs_list = ["coverage", "dop", "nsat"]

    for obj in objs:
        if obj not in objs_list:
            raise ValueError(f"Invalid objective: {obj}")

    # set tspan to 10 days
    if tspan is None:
        if dt is None:
            dt = 300.0  # [s] time step for propagation
        nstep_per_hour = int(3600 / dt)
        tspan = np.linspace(
            0, 30 * pnt.SECS_DAY, 30 * 24 * nstep_per_hour + 1
        )  # [s] time span for propagation (montly average, every 5 minutes)
    else:
        dt = tspan[1] - tspan[0]

    # set the user points
    x_user = (
        fibonacci_sphere(n_users) * pnt.R_MOON
    )  # [km] user positions on the surface of the Moon
    enu_mats = compute_enu_matrices(x_user)  # (n_users, 3, 3) PA->ENU matrices for user positions

    lat_users = np.rad2deg(
        np.arcsin(np.clip(x_user[:, 2] / pnt.R_MOON, -1, 1))
    )  # [deg] latitude of users
    users_used = []
    for phase in range(n_phase):
        users_used.append(
            np.where((lat_users >= lat_masks[phase][0]) & (lat_users <= lat_masks[phase][1]))[0]
        )

    n_var = n_walker * N_WALKER_PARAMS + N_PHASE_VEC  # number of variables
    n_const = int(n_walker * N_WALKER_CONSTRAINTS + N_GLOBAL_CONSTRAINTS + 2 * n_phase)

    # simulation flags
    compute_dop = "dop" in objs

    if compute_dop and use_wdop:
        compute_link_budget = True  # you need C/N0 for URE computation

    xl, xu = compute_variable_limits(n_walker, n_phase, n_var, float_sma)

    if dop_type_phase is None:
        dop_type_phase = ["hdop"] * n_phase  # default is HDOP for all phases

    if dop_target is None:
        dop_target = [5.0] * n_phase  # default is 5.0 for all phases

    # od error fitting coefficients
    if use_wdop:
        if float_sma:
            od_params_file = OPT_BASEPATH + "data/gridsearch_od/fit_oderr.npy"
            sma_file = OPT_BASEPATH + "data/gridsearch_od/smas.npy"
        else:
            od_params_file = OPT_BASEPATH + "data/gridsearch_od/fit_oderr_intsma.npy"
            sma_file = OPT_BASEPATH + "data/gridsearch_od/smas_ints.npy"
        ecc_file = OPT_BASEPATH + "data/gridsearch_od/eccs.npy"
        omega_file = OPT_BASEPATH + "data/gridsearch_od/omegas.npy"

        if os.path.exists(od_params_file):
            od_coeffs = np.load(od_params_file, allow_pickle=True)
            od_smas = np.load(sma_file, allow_pickle=True)
            od_eccs = np.load(ecc_file, allow_pickle=True)
            od_omegas = np.load(omega_file, allow_pickle=True)
        else:
            # run time error
            raise FileNotFoundError(
                f"OD fitting parameters not found. Please run 'python examples/ex_constellation_optimization/src/fit_od_error.py' to generate the fitting parameters."
            )

    # set config
    config = {
        "n_walker": n_walker,
        "n_phase": n_phase,
        "sat_range_phases": sat_range_phases,
        "n_var": n_var,  # walker parameters + relative phases + allocation
        "n_obj": len(objs) * n_phase,  # number of objectives
        "n_ieq_constr": n_const,  # number of inequality constraints
        "xl": xl,  # lower bounds
        "xu": xu,  # upper bounds
        "tspan": tspan,  # time span for propagation
        "x_user": x_user,  # user positions
        "users_used": users_used,
        "enu_mats": enu_mats,  # ENU matrices for user positions
        "objs": objs,  # objectives to be computed
        "elev_mask_deg": elev_mask_deg,  # elevation mask in degrees
        "cn0_thresh_user": cn0_thresh_user,  # cn0_threshold of user
        "compute_dop": compute_dop,
        "use_wdop": use_wdop,
        "dop_type_phase": dop_type_phase,  # type of DOP to be computed
        "compute_link_budget": compute_link_budget,
        "dop_target": dop_target,
        "tai0": epoch,  # initial epoch in TAI
        "max_fail_sat": max_fail_sat,
        "launch_years": launch_years,
        "eval_years": eval_years,
        "fail_model": fail_model,
        "od_coeffs": od_coeffs,  # orbit determination error fitting parameters
        "od_smas": od_smas,
        "od_eccs": od_eccs,
        "od_omegas": od_omegas,
        "float_sma": float_sma,
        "dt": dt,
    }

    return config


def fibonacci_sphere(samples=100):
    points = np.zeros((samples, 3))
    phi = np.pi * (3.0 - np.sqrt(5))  # golden angle

    for i in range(samples - 2):
        y = 1 - (i / float(samples - 2 - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points[i] = [x, y, z]

    # add south pole point
    points[samples - 1] = [0, 0, -1.0]  # south pole

    # add north pole point
    points[samples - 2] = [0, 0, 1.0]  # north pole

    return points


def print_lunanet_x(x, config):
    n_walker = config["n_walker"]
    n_const = config["n_ieq_constr"]
    n_phase = config["n_phase"]
    float_sma = config.get("float_sma", True)
    n_tot_sat = 0
    n_sats_walker = []
    sat_to_walker = []

    for i in range(n_walker):
        print("Walker {}".format(i))
        k0 = i * N_WALKER_PARAMS
        if float_sma:
            print("  sma     [km] : {}".format(x[k0]))
        else:
            print("  sma     [km] : {} (int:{})".format(int_to_sma(int(x[k0])), int(x[k0])))
        print("  ecc     []   : {}".format(x[k0 + 1]))
        print("  w       [deg]: {}".format(np.round(x[k0 + 2]) * 180 - 90))
        print("  planes  []   : {}".format(np.round(x[k0 + 3])))
        print("  nsat    []   : {}".format(np.round(x[k0 + 4])))
        print("  phase   []   : {}".format(np.round(x[k0 + 5])))
        print("  Omega0  [deg]: {}".format(np.round(x[k0 + 6] * 180 / np.pi)))
        n_sat = int(np.round(x[k0 + 3]) * np.round(x[k0 + 4]))
        n_sats_walker.append(n_sat)
        sat_to_walker.extend([i for j in range(n_sat)])

    n_sats_walker = np.array(n_sats_walker)
    n_tot_walker = np.sum(n_sats_walker)
    sat_to_walker = np.array(sat_to_walker)

    x_phase = x[n_walker * N_WALKER_PARAMS :]
    x_phase = x_phase[:n_tot_walker]
    print("x_phase:", x_phase)
    print(" ")
    for phase in range(n_phase):
        # find index
        idx = np.where(x_phase == phase)[0]

        if len(idx) == 0:
            print(f"Phase {phase} : No satellites allocated")
        else:
            # list the number of satellites for each walker
            print(f"Phase {phase}: {len(idx)} satellites")
            num_sat_walkers = {}
            for i in range(n_walker):
                num_sat_walkers[i] = 0
            for i in idx:
                if i < n_tot_walker:
                    # find a place where n_sats_walker[i] > i
                    walker = sat_to_walker[i]
                    num_sat_walkers[walker] += 1
            for i in range(n_walker):
                print("  walker {}:  {} satellites".format(i, num_sat_walkers[i]))


def int_to_sma(sma_int):
    """Convert integer sma to real sma in km"""
    T_ratios = np.array(
        [1 / 1, 2 / 3, 1 / 2, 2 / 5, 1 / 3, 2 / 7, 1 / 4]
    )  # integer multiples of 48 hr
    T_orbit = 48 * 3600 * T_ratios[sma_int]
    sma = ((T_orbit / (2 * np.pi)) ** 2 * pnt.GM_MOON) ** (1 / 3)
    return sma


def get_sim_time(x, n_walker):
    t_sim = 48 * 3600  # default 2 day
    for i in range(n_walker):
        k0 = i * N_WALKER_PARAMS
        if x[k0] < 10:  # integer sma
            sma_i = int(x[k0])
            if sma_i == 1 or sma_i == 3 or sma_i == 5:  # include 2/3, 2/5, 2/7 (2x48hr)
                t_sim = 96 * 3600
                break

    return t_sim


def create_x(walker_params, x_phase=None):
    n_walker = len(walker_params)
    n_vars = n_walker * N_WALKER_PARAMS + N_PHASE_VEC  # number of variables

    x = np.zeros(n_vars)

    n_sats = []
    for i, walker_x in enumerate(walker_params):
        # walker parameters
        # 0: semi-major axis [km] or integer index
        # 1: eccentricity []
        # 2: sign of argument of periapsis [1 or -1]
        # 3: number of planes
        # 4: number of satellites per plane
        # 5: relative phasing
        # 6: Omega0: right ascension of the ascending node [rad]
        # 7: transmittion power (P_tx)
        k0 = i * N_WALKER_PARAMS
        x[k0 : k0 + N_WALKER_PARAMS] = walker_x
        n_sats.append(x[k0 + 3] * x[k0 + 4])

    n_sats = np.array(n_sats)
    total_sats = np.sum(n_sats)

    if x_phase is None:
        # allocate each walker to each phase
        x_phase = -1 * np.ones(N_PHASE_VEC, dtype=int)
        idx0 = 0
        for i in range(n_walker):
            idxf = idx0 + int(n_sats[i])
            x_phase[idx0:idxf] = i
            idx0 = idxf
        x[n_walker * N_WALKER_PARAMS :] = x_phase
    else:
        x[n_walker * N_WALKER_PARAMS :] = x_phase

    return x


def compute_enu_matrices(x_user):
    """
    Compute the ENU matrices for the given user positions.

    Parameters:
    x_user : np.ndarray [n_users, 3]
        Cartesian coordinates of user positions.
    Returns:
    enu_mats : np.ndarray [n_users, 3, 3]
        ENU matrices for each user position.
    """
    enu_mats = np.zeros((x_user.shape[0], 3, 3))

    phi = np.arcsin(np.clip(x_user[:, 2] / pnt.R_MOON, -1, 1))  # [rad] latitudes
    lam = np.arctan2(x_user[:, 1], x_user[:, 0])  # [rad] longitudes

    enu_mats[:, 0, 0] = -np.sin(lam)  # East
    enu_mats[:, 0, 1] = np.cos(lam)
    enu_mats[:, 0, 2] = 0  # North
    enu_mats[:, 1, 0] = -np.sin(phi) * np.cos(lam)  # North
    enu_mats[:, 1, 1] = -np.sin(phi) * np.sin(lam)
    enu_mats[:, 1, 2] = np.cos(phi)  # Up
    enu_mats[:, 2, 0] = np.cos(phi) * np.cos(lam)  # Up
    enu_mats[:, 2, 1] = np.cos(phi) * np.sin(lam)
    enu_mats[:, 2, 2] = np.sin(phi)  # East

    return enu_mats


def compute_variable_limits(n_walker, n_phase, n_var, float_sma):

    # walker parameters
    # 0: semi-major axis [km]
    # 1: eccentricity []
    # 2: sign of argument of periapsis [1 or -1]
    # 3: number of planes
    # 4: number of satellites per plane
    # 5: relative phase
    # 6: Omega0: right ascension of the ascending node [rad]
    # 7: transmittion power (P_tx)
    xl = np.zeros(n_var)
    xu = np.inf * np.ones(n_var)

    k = 0  # index
    for i in range(n_walker):
        # 0 semi-major axis
        # xl[k] = 3474
        # xu[k] = 17370
        if float_sma:
            xl[k] = 3500
            xu[k] = 13000
        else:
            xl[k] = 0
            xu[k] = N_SMA_INT - 1
        k += 1
        # 1 eccentricity
        xl[k] = 0
        xu[k] = 0.7
        k += 1
        # 2 sign of argument of periapsis
        xl[k] = 0
        xu[k] = 1
        k += 1
        # 3 number of planes
        xl[k] = 2
        xu[k] = 6
        k += 1
        # 4 number of satellites per plane
        xl[k] = 2
        xu[k] = 8
        k += 1
        # 5 relative phase
        xl[k] = 1
        xu[k] = 2
        k += 1
        # 6 Omega 0
        xl[k] = 0.0
        xu[k] = 2 * np.pi
        k += 1

    # phase allocation variable
    for i in range(k, n_var):
        xl[i] = -1
        xu[i] = n_phase - 1

    return xl, xu


def compute_constraints(x, config):
    """Compute the constraints for the constellation optimization problem. g(x) <= 0

    Parameters
    ----------
    x : _type_
        _description_
    config : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    n_walker = config["n_walker"]
    n_const = config["n_ieq_constr"]
    n_phase = config["n_phase"]
    sat_phase_ranges = config["sat_range_phases"]
    float_sma = config.get("float_sma", True)

    consts = np.zeros(n_const)
    j = 0
    n_sats_walker = 0

    # constsnts
    R_MOON = 1737.4
    MIN_ALT = 100.0

    for i in range(n_walker):
        k0 = i * N_WALKER_PARAMS
        # constraint: periapsis altitude constraint > 0
        if float_sma:
            a = x[k0]
        else:
            a = int_to_sma(int(x[k0]))
        e = x[k0 + 1]
        consts[j] = -(a * (1 - e) - R_MOON - MIN_ALT)
        j += 1

        # constraint: phasing < number of planes > 0
        phasing = x[k0 + 5]
        planes = x[k0 + 3]
        consts[j] = phasing - planes
        j += 1

        # add nsats
        n_sat_plane = x[k0 + 4]
        n_sats_walker += planes * n_sat_plane

    # number of used satellites
    x_phase = x[n_walker * N_WALKER_PARAMS :]  # phase allocation
    x_phase = x_phase[: int(n_sats_walker)]  # only consider allocated satellites
    n_sats_used = np.sum(x_phase >= 0)
    n_sat_init = np.sum(x_phase == 0)  # number of initial satellites

    # constraint: maximum number of satellites should not exceed allower number
    k = j
    consts[k] = n_sats_used - n_sats_walker
    k += 1
    consts[k] = n_sats_walker - N_PHASE_VEC
    k += 1

    n_sat_phase = 0  # number of satellites in each phase
    for i in range(n_phase):
        n_sat_phase += np.sum(x_phase == i)
        consts[k] = sat_phase_ranges[i][0] - n_sat_phase
        k += 1
        consts[k] = n_sat_phase - sat_phase_ranges[i][1]
        k += 1
        # print("Phase {}:  n_sat = {}   (min: {}, max: {})".format(i, n_sat_phase, sat_phase_ranges[i][0], sat_phase_ranges[i][1]))

    params = {}
    params["n_sat_used"] = n_sats_used
    params["n_sats_walker"] = n_sats_walker
    params["n_sat_init"] = n_sat_init
    params["n_sat_walker"] = n_sats_walker

    return np.array(consts), params


def setup_hybrid_walker_constellation(
    x, et, dyn, float_sma, n_walker=1, n_phase=3, compute_link_budget=False, debug=False
):
    """
    Set up Hybrid Walker Constellation based on the provided parameters.
    """

    out = {}  # output parameters

    # walker parameters
    # 0: semi-major axis [km]
    # 1: eccemtrocototu
    # 2: sign of argument of periapsis [1 or -1]
    # 3: number of planes
    # 4: number of satellites per plane
    # 5: relative phase
    # 6: Omega0: right ascension of the ascending node [rad]
    # 7: transmittion power (P_tx)
    n_param_walker = N_WALKER_PARAMS  # number of parameters for each walker
    nx_walker = n_walker * n_param_walker

    x_walker = x[:nx_walker]
    x_alloc = x[nx_walker:]  # allocation of satellites to phases

    coes = np.zeros((0, 6))  # Initialize empty array for orbital elements
    walker_idxs = []
    lunanet_antennas = []
    for i in range(n_walker):
        x_walker_i = x_walker[i * n_param_walker : i * n_param_walker + n_param_walker]
        coes_tmp = setup_walker(x_walker_i, float_sma, debug=debug)
        P_tx_w = P_TX_W
        coes = np.vstack((coes, coes_tmp))
        n_sat_walker = coes_tmp.shape[0]
        if n_sat_walker > 0:
            walker_idxs.extend([i for j in range(n_sat_walker)])
            if compute_link_budget:
                lunanet_antenna_walker = LunanetSatAntenna(P_tx_w, coes_tmp[0])
                if debug:
                    print("Lunanet Antenna for Walker {}: P_tx = {}".format(i, P_tx_w))
                    lunanet_antenna_walker.plot_pattern()
                lunanet_antennas.extend([lunanet_antenna_walker for i in range(n_sat_walker)])

    walker_idxs = np.array(walker_idxs)
    n_sat_walker = coes.shape[0]
    n_sat_used = np.sum(x_alloc >= 0)  # number of satellites used in the constellation

    x_alloc = x_alloc[:n_sat_walker]  # filter out redundant satellites

    # propagate orbits
    coes = coes[x_alloc >= 0, :]  # filter out unused satellites

    x_orb_pa, x_orb_mci = propagate_lunanet_orbits(coes, et, dyn)

    if compute_link_budget:
        lunanet_antennas = [lunanet_antennas[i] for i in range(n_sat_walker) if x_alloc[i] >= 0]

    x_phase = x_alloc[x_alloc >= 0]  # filter out unused satellites
    walker_idxs = walker_idxs[x_alloc >= 0]  # filter out unused satellites

    out["x_orb_mci"] = x_orb_mci
    out["n_sat_used"] = n_sat_used
    out["coes"] = coes
    out["walker_idxs"] = walker_idxs

    return x_orb_pa, x_phase, lunanet_antennas, out


def propagate_lunanet_orbits(coes, et, dyn):
    """
    Propagate the orbits using two-body dynamics.

    Parameters:
    coes : np.ndarray
        Classical orbital elements of shape (n_sat, 6).
    tspan : np.ndarray
        Time span for propagation.
    dt : float
        Time step for propagation.

    Returns:
    x_orb : np.ndarray
        Propagated orbits in Cartesian coordinates.
    """
    n_sat = coes.shape[0]
    lent = len(et)
    x_orb_pa = np.zeros((n_sat, lent, 6))
    x_orb_ci = np.zeros((n_sat, lent, 6))

    for si in range(n_sat):
        rv0_op = pnt.classical_to_cart(coes[si, :], pnt.GM_MOON)
        rv0_mci = pnt.convert_frame(et[0], rv0_op, pnt.MOON_OP, pnt.MOON_CI)
        x_orb_ci[si] = dyn.propagate(rv0_mci, et[0], et)
        x_orb_pa[si] = pnt.convert_frame(et, x_orb_ci[si], pnt.MOON_CI, pnt.MOON_PA)

    return x_orb_pa, x_orb_ci


def setup_walker(x_walker, float_sma, debug=False):
    """
    Set up a Walker Constellation based on the provided parameters.
    """
    if float_sma:
        sma = x_walker[0]  # semi-major axis [km]
    else:
        sma = int_to_sma(int(x_walker[0]))  # semi-major axis [km]

    ecc = x_walker[1]  # eccentricity []
    wsign = np.round(x_walker[2])  # argument of periapsis [rad]
    p = np.round(x_walker[3])  # number of planes [-]
    t = np.round(x_walker[4])  # number of satellites per plane [-]
    f = np.round(x_walker[5])  # relative phase [-]
    Omega0 = x_walker[6]  # right ascension of the ascending node [rad]

    if wsign > 0.5:
        w = 90 * pnt.RAD  # argument of periapsis [rad]
    else:
        w = -90 * pnt.RAD

    cosinc = np.sqrt((1 - ecc**2) * 3 / 5)  # inclination for frozen orbit
    inc = np.arccos(cosinc)

    # Initialize list to store orbits
    tot_sat = t * p
    coes = np.zeros((int(t * p), 6))

    if debug:
        print("Walker Constellation Parameters:")
        print("  sma        [km] : {}".format(sma))
        print("  ecc        []   : {}".format(ecc))
        print("  inc        [deg]: {}".format(np.round(inc * 180 / np.pi)))
        print("  w          [deg]: {}".format(np.round(w * 180 / np.pi)))
        print("  planes     []   : {}".format(int(p)))
        print("  sat/plane  []   : {}".format(int(t)))
        print("  rel phase  []   : {}".format(int(f)))
        print("  Omega0     [deg]: {}".format(np.round(Omega0 * 180 / np.pi)))
        print("  tot_sat    []   : {}".format(int(tot_sat)))

    idx = 0
    for plane in range(int(p)):
        for s in range(int(t)):
            alpha = s * 2 * np.pi / t  # intra-plane spacing
            beta = plane * f * 2 * np.pi / tot_sat  # inter-plane spacing
            M = pnt.wrap2pi(alpha + beta)  # true anomaly
            Omega = pnt.wrap2pi(
                Omega0 + plane * 2 * np.pi / p
            )  # right ascension of the ascending node
            # M = pnt.true_to_mean_anomaly(nu, ecc)  # mean anomaly
            coes[idx] = np.array([sma, ecc, inc, Omega, w, M])
            idx += 1

    return coes


def compute_objectives(x, config, debug=False):
    """
    Compute the objectives for the constellation optimization problem.
    """

    params = {}  # output parameter

    if debug:
        time_start = time.time()

    sma_float = config.get("float_sma", True)  # optimize semi-major axis
    if not sma_float:
        tspan = config["tspan"]
        lent = tspan.shape[0]
    else:
        dt = config.get("dt", config["tspan"][1] - config["tspan"][0])
        tsim = get_sim_time(x, config["n_walker"])
        tspan = np.linspace(0, tsim, int(tsim / dt) + 1)
        lent = tspan.shape[0]

    if debug:
        print("  lent:", lent)

    n_walker = config["n_walker"]
    n_phase = config["n_phase"]
    et0 = config["tai0"]
    objs = config["objs"]
    dop_target = config["dop_target"]
    float_sma = config.get("float_sma", True)
    compute_link_budget = config.get("compute_link_budget", False)

    # dynamics
    dyn = pnt.NBodyDynamics()
    dyn.set_integrator(pnt.IntegratorType.RKF45)
    dyn.set_integrator_params(pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10))
    dyn.add_body(pnt.Body.Moon(2, 2))
    dyn.add_body(pnt.Body.Earth())
    dyn.add_body(pnt.Body.Sun())
    dyn.set_time_step(60)
    dyn.set_frame(pnt.MOON_CI)

    t_tai = et0 + tspan  # [s] time since first epoch in TAI

    if debug:
        print("  Time to set up dynamics: {:.2f} s".format(time.time() - time_start))
        time_start = time.time()

    # Setup Hybrid Walker Constellation
    x_orb, x_phase, lunanet_antennas, params_walker = setup_hybrid_walker_constellation(
        x,
        et0 + tspan,
        dyn,
        float_sma,
        n_walker,
        n_phase,
        compute_link_budget=compute_link_budget,
        debug=False,
    )
    nsat_used = x_orb.shape[0]
    params["n_sat_used"] = nsat_used
    if debug:
        print("  Time to set up constellation: {:.2f} s".format(time.time() - time_start))
        time_start = time.time()

    # Compute coverage and DOP
    coverage, dop, ures_noise, od_err, probs = compute_coverage_dop(
        t_tai,
        x_orb,
        x_phase,
        lunanet_antennas=lunanet_antennas,
        coe=params_walker["coes"],
        config=config,
    )
    if debug:
        print("  Time to compute coverage and DOP: {:.2f} s".format(time.time() - time_start))
        time_start = time.time()

    # Compute objectives based on the requested ones
    results = []
    result_labels = []
    for phase in range(n_phase):
        for obj in objs:
            if obj == "dop":
                dop_phase = dop[phase]  # (n_user, lent, case)
                n_user = dop_phase.shape[0]
                prob_phase = probs[phase]
                dop_target_ratio = (
                    np.sum(np.sum(dop_phase <= dop_target[phase], axis=0), axis=0) / lent / n_user
                )  # (n_user, case)
                wsum_dop_target_ratio = np.sum(prob_phase * dop_target_ratio)
                objval = 1 - wsum_dop_target_ratio
                results.append(objval)
                result_labels.append(f"pdop_p{phase}")
            elif obj == "nsat":
                nsat = np.sum(x_phase <= phase)
                nsat_max = config["sat_range_phases"][phase][1]
                nsat_min = config["sat_range_phases"][phase][0]
                objval = (nsat - nsat_min) / (nsat_max - nsat_min)
                results.append(objval)
                result_labels.append(f"nsat_p{phase}")
            else:
                raise ValueError(f"Unknown objective: {obj}")

    return np.array(results), result_labels, params
