from src.navigation import compute_dop, compute_gnss_pseudorange_noise, DOP_INF
from src.failure_model import NavSatFailureModel, compute_fail_probs
from src.communication import LunanetSatAntenna, compute_lunanet_cn0
from src.navigation import compute_gnss_pseudorange_noise
from src.odsim import param_to_oderr
import numpy as np
import pylupnt as pnt
from src.navigation import sim_covdop_nfail


def compute_coverage_dop(t_tai, x_orb, x_phase, lunanet_antennas, config, coe, debug=False):
    """
    Compute the coverage and DOP for a given constellation and user positions.

    Returns:
    coverage_phases : list of np.ndarray [n_user, lent]
    dop_phases : list of np.ndarray [n_user, lent]
    """

    min_elev_deg = config["elev_mask_deg"]
    cn0_thresh = config["cn0_thresh_user"]
    use_weighted_dop = config["use_wdop"]
    enu_mats = config["enu_mats"]  # (n_user, 3, 3)
    compute_link_budget = config["compute_link_budget"]
    x_user = config["x_user"]  # (n_user, 3)
    n_phase = config["n_phase"]
    dop_target = config["dop_target"]

    n_sat, lent, _ = x_orb.shape
    n_user = x_user.shape[0]

    # Compute user-to-satellite vectors
    r_m2u = np.tile(x_user[:, None, None, :], (1, n_sat, lent, 1))  # (n_user, n_sat, lent, 3)
    r_m2s = np.tile(x_orb[None, :, :, :3], (n_user, 1, 1, 1))  # (n_user, n_sat, lent, 3)
    r_u2s = r_m2s - r_m2u  # (n_user, n_sat, lent, 3)

    # Convert to ENU
    r_u2s_enu = np.matmul(enu_mats[:, None, None, :, :], r_u2s[..., None]).squeeze(-1)
    r_u2s_v = r_u2s_enu.reshape(-1, 3)
    u2s_norm = np.linalg.norm(r_u2s_v, axis=1)
    r_u2s_v /= u2s_norm[:, None]
    r_u2s_enu = r_u2s_v.reshape(n_user, n_sat, lent, 3)

    elev_user = np.arcsin(r_u2s_enu[..., 2])  # elevation = arcsin(z)
    coverage = elev_user >= np.deg2rad(min_elev_deg)  # (n_user, n_sat, lent)

    # Compute C/N0 -------------------------------------------------------------
    if compute_link_budget:
        if debug:
            print("Compute C/N0...")
        cn0_links, coverage = compute_lunanet_cn0(
            t_tai,
            x_orb,
            r_u2s,
            elev_user,
            u2s_norm.reshape(n_user, n_sat, lent),
            lunanet_antennas,
            min_elev_deg,
            cn0_thresh,
        )
    else:
        cn0_links = cn0_thresh * np.ones((n_user, n_sat, lent))

    # Compute URE -------------------------------------------------------------
    ures_noise = np.ones_like(cn0_links)
    if use_weighted_dop:
        if debug:
            print("Compute URE based on C/N0...")
        ures_noise = compute_gnss_pseudorange_noise(cn0_links)  # (n_user, n_sat, lent)

    # Do OD Simulations -------------------------------------------------------------
    if use_weighted_dop:
        od_errs = get_oderrs_from_log(coe, t_tai, config, debug)  # (n_sat, lent)
    else:
        od_errs = np.zeros((n_sat, lent))

    # add od error to ures
    ure_others = 3.0  # other sources of URE [m]
    ures = ures_noise + od_errs[None, :, :] + ure_others  # (n_user, n_sat, lent)

    # DOP and coverage for each phase ------------------------------------------------
    coverage_phases = []
    dop_phases = []
    prob_phases = []

    for phase in range(n_phase):
        # Failure probabilities
        fail_probs = compute_fail_probs(phase, x_phase, n_sat, config)
        # fail_prob_sat = 0.01  # placeholder for failure probability of each satellite
        # fail_probs = np.full(n_sat, fail_prob_sat)  # placeholder

        # simulation config
        max_fail = config[
            "max_fail_sat"
        ]  # maximum number of failures allowed in a visibility pattern

        if debug:
            print(f"Simulate coverage and DOP for phase {phase} with max {max_fail} failures...")
        cov_phase, dop_phase, probs_phase = sim_covdop_nfail(
            r_u2s_enu, ures, coverage, x_phase, max_fail, fail_probs, phase, config, debug=debug
        )

        if debug:
            print("Phase: ", phase + 1)
            n_cases = dop_phase.shape[2] - 1
            n_user_phase = dop_phase.shape[0]
            print("  Failure cases: ", n_cases)
            dop_target_ratio10 = (
                np.sum(np.sum(dop_phase <= 10, axis=0), axis=0) / lent / n_user_phase
            )  # (case)
            dop_target_ratio15 = (
                np.sum(np.sum(dop_phase <= 15, axis=0), axis=0) / lent / n_user_phase
            )  # (case)
            dop_target_ratio20 = (
                np.sum(np.sum(dop_phase <= 20, axis=0), axis=0) / lent / n_user_phase
            )  # (case)
            print(
                "  DOP no error:  (<10) : {}  (<15): {}  (<20) : {}".format(
                    dop_target_ratio10[0], dop_target_ratio15[0], dop_target_ratio20[0]
                )
            )
            for i in range(1, n_cases + 1):
                print(
                    "    case {}: prob={:.4f}  (<10) : {:.4f}  (<15): {:.4f}  (<20) : {:.4f}".format(
                        i,
                        probs_phase[i],
                        dop_target_ratio10[i],
                        dop_target_ratio15[i],
                        dop_target_ratio20[i],
                    )
                )
            print(" ")

        coverage_phases.append(cov_phase)
        dop_phases.append(dop_phase)
        prob_phases.append(probs_phase)

    return coverage_phases, dop_phases, ures_noise, od_errs, prob_phases


def get_oderrs_from_log(coe, t_tai, config, debug):
    if debug:
        print("Compute OD errors...")

    n_sat = coe.shape[0]
    lent = t_tai.shape[0]

    T_orbits = 2 * np.pi * np.sqrt(coe[:, 0] ** 3 / pnt.GM_MOON)  # [s] orbital periods
    od_errs = np.zeros((n_sat, lent))  # [m] orbit determination errors

    for i in range(n_sat):
        t_in_period = (t_tai - t_tai[0]) % T_orbits[i]  # time within the orbital period
        M = (2 * np.pi / T_orbits[i]) * t_in_period + coe[i, 5]  # mean anomaly
        # brind M to [0, 2pi]
        M = pnt.wrap_to_two_pi(M)

        sma = coe[i, 0]
        ecc = coe[i, 1]
        Omega = coe[i, 3]

        # find the closest sma, ecc, omega in the od grid
        idx_sma = np.argmin(np.abs(config["od_smas"] - sma))
        idx_ecc = np.argmin(np.abs(config["od_eccs"] - ecc))
        idx_omega = np.argmin(np.abs(config["od_omegas"] - Omega))
        params = config["od_coeffs"][idx_sma, idx_ecc, idx_omega]

        od_err = param_to_oderr(M, params)  # [m] orbit determination error
        od_errs[i, :] = od_err * 1000  # convert to meters

    return od_errs
