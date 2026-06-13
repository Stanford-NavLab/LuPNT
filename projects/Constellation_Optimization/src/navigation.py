import numpy as np
from numba import njit
from src.communication import LunanetReceiverParam

import numpy as np
from numba import njit, prange

# ---------------- DOP helpers (from previous message) ----------------
DOP_INF = 100.0


@njit(cache=True, fastmath=True)
def _compute_pdop_core(H, w_diag):
    n_vis = H.shape[0]
    if n_vis < 4:
        return DOP_INF
    Hw = H * w_diag.reshape(n_vis, 1)
    M = Hw.T @ H
    eps = 1e-12
    for i in range(4):
        M[i, i] += eps
    Minv = np.linalg.inv(M)
    tr3 = Minv[0, 0] + Minv[1, 1] + Minv[2, 2]
    if tr3 <= 0.0:
        return DOP_INF
    return np.sqrt(tr3)


@njit(cache=True, fastmath=True)
def _compute_hdop_dem_core(H, w_diag, dem_sigma):
    n_vis = H.shape[0]
    if n_vis < 3:
        return DOP_INF
    Hw = H * w_diag.reshape(n_vis, 1)
    M = Hw.T @ H
    # add DEM pseudo-measurement directly to z
    w_dem = 1.0 / (dem_sigma * dem_sigma)
    M[2, 2] += w_dem
    eps = 1e-12
    for i in range(4):
        M[i, i] += eps
    Minv = np.linalg.inv(M)
    s = Minv[0, 0] + Minv[1, 1]
    if s <= 0.0:
        return DOP_INF
    return np.sqrt(s)


# ---------------- Core Monte Carlo (JIT) ----------------
@njit(cache=True, fastmath=True, parallel=True)
def _sim_covdop_nfail_core(
    r_u2s_enu,
    ures,
    coverage,
    used_mask,
    sat_used_idx,
    max_fail,
    dop_type_flag,
    use_wdop_flag,
    compute_dop_flag,
    dem_sigma,
):
    """
    r_u2s_enu : (n_user, n_sat_all, lent, 3) float64
    ures      : (n_user, n_sat_all, lent) float64
    coverage  : (n_user, n_sat_all, lent) bool
    used_mask : (n_sat_all,) bool  (x_phase <= phase)
    sat_used_idx : (n_sat,) int64
    dop_type_flag : 0=PDOP, 1=HDOP+DEM
    """
    n_user, n_sat_all, lent, _ = r_u2s_enu.shape
    n_sat = sat_used_idx.shape[0]

    # number of scenarios (no fail, + single-fail, + double-fail)
    n_sim = 1
    if max_fail >= 1:
        n_sim += n_sat
    if max_fail >= 2:
        n_sim += n_sat * (n_sat - 1) // 2

    cov_phase = np.zeros((n_user, lent, n_sim), dtype=np.float64)
    dop_phase = np.empty((n_user, lent, n_sim), dtype=np.float64)

    min_sats = 4 if dop_type_flag == 0 else 3

    for i in prange(n_user):
        for t in range(lent):

            # base visibility (no failures)
            vis0 = np.zeros(n_sat_all, dtype=np.uint8)
            n_vis0 = 0
            for s in range(n_sat_all):
                if used_mask[s] and coverage[i, s, t]:
                    vis0[s] = 1
                    n_vis0 += 1

            # If even the no-fail case is under-minimum, all fail cases are too
            if n_vis0 < min_sats or compute_dop_flag == 0:
                for m in range(n_sim):
                    cov_phase[i, t, m] = 0.0
                    dop_phase[i, t, m] = DOP_INF
                continue

            # buffers reused across scenarios
            Hbuf = np.empty((n_sat_all, 4), dtype=np.float64)
            wbuf = np.empty(n_sat_all, dtype=np.float64)

            # helper to build H/w for a given failed pair (a,b), where -1=none
            def build_and_eval(fail_a, fail_b):
                row = 0
                for s in range(n_sat_all):
                    if vis0[s] == 0:
                        continue
                    if s == fail_a or s == fail_b:
                        continue
                    r = r_u2s_enu[i, s, t]
                    Hbuf[row, 0] = r[0]
                    Hbuf[row, 1] = r[1]
                    Hbuf[row, 2] = r[2]
                    Hbuf[row, 3] = 1.0
                    if use_wdop_flag == 1:
                        sig = ures[i, s, t]
                        if sig < 1e-12:
                            sig = 1e-12
                        wbuf[row] = 1.0 / (sig * sig)
                    else:
                        wbuf[row] = 1.0
                    row += 1

                if row < min_sats:
                    return 0.0, DOP_INF

                if dop_type_flag == 0:
                    dop_val = _compute_pdop_core(Hbuf[:row, :], wbuf[:row])
                else:
                    dop_val = _compute_hdop_dem_core(Hbuf[:row, :], wbuf[:row], dem_sigma)
                return 1.0, dop_val

            # 0) no-failure
            cov, dop = build_and_eval(-1, -1)
            cov_phase[i, t, 0] = cov
            dop_phase[i, t, 0] = dop

            # 1) single failures
            if max_fail >= 1:
                for k in range(n_sat):
                    fail_a = sat_used_idx[k]
                    cov, dop = build_and_eval(fail_a, -1)
                    idx = k + 1
                    cov_phase[i, t, idx] = cov
                    dop_phase[i, t, idx] = dop

            # 2) double failures
            if max_fail >= 2:
                m = n_sat + 1
                for k in range(n_sat):
                    a = sat_used_idx[k]
                    for l in range(k + 1, n_sat):
                        b = sat_used_idx[l]
                        cov, dop = build_and_eval(a, b)
                        cov_phase[i, t, m] = cov
                        dop_phase[i, t, m] = dop
                        m += 1

    return cov_phase, dop_phase


# ---------------- Thin Python wrapper (dicts/strings handled here) ----------------
def sim_covdop_nfail(
    r_u2s_enu,
    ures,
    coverage,
    x_phase,
    max_fail,
    fail_probs,
    phase,
    config,
    debug=False,
    dem_sigma=1.0,
):
    """
    Wrapper that extracts config & masks, computes scenario probabilities once in Python,
    and calls the Numba core for the heavy lifting.
    """
    users_used_idx = np.asarray(config["users_used"][phase], dtype=np.int64)
    dop_type_flag = 0 if config["dop_type_phase"][phase] == "pdop" else 1
    use_wdop_flag = 1 if config["use_wdop"] else 0
    compute_dop_flag = 1 if config["compute_dop"] else 0

    # Subset users up front (Numba doesn't like fancy dicts/logic)
    r_u2s_enu_u = np.ascontiguousarray(r_u2s_enu[users_used_idx, :, :, :], dtype=np.float64)
    ures_u = np.ascontiguousarray(ures[users_used_idx, :, :], dtype=np.float64)
    coverage_u = np.ascontiguousarray(coverage[users_used_idx, :, :], dtype=np.bool_)

    used_mask = x_phase <= phase  # (n_sat_all,) bool
    sat_used_idx = np.where(used_mask)[0].astype(np.int64)

    # Run the JIT core
    cov_phase, dop_phase = _sim_covdop_nfail_core(
        r_u2s_enu_u,
        ures_u,
        coverage_u,
        used_mask.astype(np.bool_),
        sat_used_idx,
        int(max_fail),
        int(dop_type_flag),
        int(use_wdop_flag),
        int(compute_dop_flag),
        float(dem_sigma),
    )

    # ---- scenario probabilities (Python, done once) ----
    # Build probs for used satellites only, in the same order as sat_used_idx
    fp = np.asarray(fail_probs, dtype=np.float64)[sat_used_idx]
    n_sat = fp.size

    n_sim = 1
    if max_fail >= 1:
        n_sim += n_sat
    if max_fail >= 2:
        n_sim += n_sat * (n_sat - 1) // 2

    probs = np.zeros(n_sim, dtype=np.float64)

    # 0) no-failure
    p0 = 1.0
    for j in range(n_sat):
        p0 *= 1.0 - fp[j]
    probs[0] = p0

    # 1) single-failure
    idx = 1
    if max_fail >= 1:
        for k in range(n_sat):
            pk = fp[k]
            for j in range(n_sat):
                if j != k:
                    pk *= 1.0 - fp[j]
            probs[idx] = pk
            idx += 1

    # 2) double-failure
    if max_fail >= 2:
        for k in range(n_sat):
            for l in range(k + 1, n_sat):
                pkl = fp[k] * fp[l]
                for j in range(n_sat):
                    if j != k and j != l:
                        pkl *= 1.0 - fp[j]
                probs[idx] = pkl
                idx += 1

    total_prob = probs.sum()
    probs_phase = probs / total_prob if total_prob > 0.0 else probs

    return cov_phase, dop_phase, probs_phase


# ---------- public API (python dispatch) ----------
def compute_dop(Hmat, ures, dop_type, use_wdop=False, dem_sigma=1.0):
    """
    Hmat : (n_visible, 4) geometry matrix
    ures : (n_visible,) user range errors (σ_i) in meters
    dop_type : "pdop" or "hdop"
    use_wdop : if True, weight by 1/σ_i^2; else all ones
    dem_sigma : std dev [m] for DEM constraint (used only for "hdop")
    """
    DOP_INF = 100.0
    H = np.ascontiguousarray(Hmat, dtype=np.float64)

    if use_wdop:
        # w_i = 1 / σ_i^2 ; guard zeros
        ures = np.asarray(ures, dtype=np.float64)
        w_diag = 1.0 / np.maximum(ures, 1e-12) ** 2
    else:
        w_diag = np.ones(H.shape[0], dtype=np.float64)

    if dop_type == "pdop":
        return _compute_pdop_core(H, w_diag)
    elif dop_type == "hdop":
        return _compute_hdop_dem_core(H, w_diag, float(dem_sigma))
    else:
        return DOP_INF


def compute_gnss_pseudorange_noise(CN0_dB):

    receiver_param = LunanetReceiverParam()

    CN0 = 10 ** (CN0_dB / 10)
    Bn = receiver_param.Bn
    Rc = receiver_param.Rc  # Chip rate in Hz
    Bfe = receiver_param.b * Rc
    T = receiver_param.T
    D = receiver_param.D
    Tc = 1.0 / Rc
    C = 299792458  # Speed of light in m/s

    sigma = np.zeros_like(CN0)

    case1 = D >= (np.pi * Rc / Bfe)
    case2 = (D > (Rc / Bfe)) & (~case1)
    case3 = ~case1 & ~case2

    # Case 1: Wide spacing discriminator
    if case1:
        sigma = np.sqrt(Bn / (2.0 * CN0) * D * (1.0 + 2.0 / (T * CN0 * (2 - D))))
    if case2:
        tmp1 = Bn / (2.0 * CN0)
        tmp2 = 1.0 / (Bfe * Tc) + Bfe * Tc / (np.pi - 1) * (D - 1.0 / (Bfe * Tc)) ** 2
        tmp3 = 1.0 + 2.0 / (T * CN0 * (2 - D))
        sigma = np.sqrt(tmp1 * tmp2 * tmp3)
    if case3:
        sigma = np.sqrt(Bn / (2.0 * CN0) * (1.0 / (Bfe * Tc)) * (1.0 + 1.0 / (T * CN0)))

    sigma_m = sigma * (C * Tc)

    return sigma_m
