import numpy as np
import pylupnt as pnt
from tqdm import tqdm
from src.communication import LunanetSatAntenna, compute_lunanet_cn0
from src.navigation import compute_gnss_pseudorange_noise
from concurrent.futures import ProcessPoolExecutor, as_completed
from src.constellation_design import fibonacci_sphere, compute_enu_matrices, setup_walker
import os
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pickle
from matplotlib.colors import LogNorm
from numba import njit, prange
import math
import matplotlib.patheffects as pe


def run_one_coverage_case(args):
    (
        i,
        j,  # indices
        sma,
        inc,
        walker_pattern,  # parameters
        et0,
        dt_sim,
        min_elev_deg,
        sim_num,
        sim_id,
        x_user,
        enu_mats,
        pole_mask,
        hybrid,
    ) = args

    import pylupnt as pnt  # import inside worker to avoid pickling module

    # local helpers assumed importable; if they live in this file, that's fine
    # otherwise, import from your package:
    # from yourpkg.whatever import fibonacci_sphere, compute_enu_matrices, setup_walker

    pnt.set_lupnt_epoch(0.0)

    n_res = 2

    n_planes, n_planes_sat = walker_pattern

    # Skip infeasible eccentricities
    ecc2 = 1 - 5 / 3 * np.cos(inc) ** 2
    if ecc2 < 0:
        return (i, j, np.nan, np.nan)

    ecc = float(np.sqrt(max(0.0, ecc2)))
    peri_h = sma * (1 - ecc) - pnt.R_MOON
    if peri_h < 100.0:
        return (i, j, np.nan, np.nan)

    # deploy satellites
    wsign = 1  # w = 90 deg
    f = 1  # f = 1
    Omega0 = 0.0

    if not hybrid:
        x_walker = np.array([sma, ecc, wsign, n_planes, n_planes_sat, f, Omega0], dtype=float)
        coes = setup_walker(x_walker, float_sma=True)
    else:
        # combine north and south walker constellations
        wsign = 1
        x_walker_s = np.array([sma, ecc, wsign, n_planes, n_planes_sat, f, Omega0], dtype=float)
        coes_s = setup_walker(x_walker_s, float_sma=True)

        wsign = 0
        x_walker_n = np.array([sma, ecc, wsign, n_planes, n_planes_sat, f, Omega0], dtype=float)
        coes_n = setup_walker(x_walker_n, float_sma=True)

        coes = np.vstack((coes_s, coes_n))

    n_sat = coes.shape[0]

    # dynamics
    dynamics = pnt.NBodyDynamics()
    dynamics.set_integrator(pnt.IntegratorType.RKF45)
    dynamics.set_integrator_params(pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10))
    dynamics.add_body(pnt.Body.Moon(2, 2))
    dynamics.add_body(pnt.Body.Earth())
    dynamics.add_body(pnt.Body.Sun())
    dynamics.set_time_step(dt_sim)
    dynamics.set_frame(pnt.MOON_CI)

    # one orbit period
    T_sim = 27.321661 * 24 * 3600  # [s] one sidereal month
    T_orbit = 2 * np.pi * np.sqrt(sma**3 / pnt.GM_MOON)
    n_orbits = max(1, int(T_sim / T_orbit) + 1)  # simulate at least one full orbit
    n_sim = int(T_orbit * n_orbits / dt_sim)
    tspan = np.linspace(0.0, T_orbit * n_orbits, n_sim + 1)
    et = et0 + tspan
    lent = len(tspan)

    # propagate orbits (Moon-CI → Moon-PA)
    x_orb = np.zeros((n_sat, lent, 6), dtype=float)
    for si in range(n_sat):
        rv0_op = pnt.classical_to_cart(coes[si, :], pnt.GM_MOON)
        rv0_mci = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)
        x_sat_mci = dynamics.propagate(rv0_mci, et)
        x_orb[si] = pnt.convert_frame(et, x_sat_mci, pnt.MOON_CI, pnt.MOON_PA)

    # user-to-satellite vectors
    n_user = x_user.shape[0]
    r_m2u = np.tile(x_user[:, None, None, :], (1, n_sat, lent, 1))  # (n_user, n_sat, lent, 3)
    r_m2s = np.tile(x_orb[None, :, :, :3], (n_user, 1, 1, 1))  # (n_user, n_sat, lent, 3)
    r_u2s = r_m2s - r_m2u

    # ENU conversion
    r_u2s_enu = np.matmul(enu_mats[:, None, None, :, :], r_u2s[..., None]).squeeze(-1)
    r_u2s_v = r_u2s_enu.reshape(-1, 3)
    u2s_norm = np.linalg.norm(r_u2s_v, axis=1)
    # avoid divide-by-zero if any exact overlaps occur
    u2s_norm[u2s_norm == 0] = 1.0
    r_u2s_v /= u2s_norm[:, None]
    r_u2s_enu = r_u2s_v.reshape(n_user, n_sat, lent, 3)

    elev_user = np.arcsin(r_u2s_enu[..., 2])
    coverage = elev_user >= np.deg2rad(min_elev_deg)
    nsat_obs = np.sum(coverage, axis=1)  # (n_user, lent)

    coverage_4fold = np.sum(nsat_obs >= 4) / (n_user * lent)

    # pole coverage
    n_pole = int(np.sum(pole_mask))
    if n_pole == 0:
        pole_coverage_4fold = np.nan
    else:
        nsat_obs_pole = nsat_obs[pole_mask, :]
        pole_coverage_4fold = np.sum(nsat_obs_pole >= 4) / (n_pole * lent)

    print(
        f"Completed case {sim_id}/{sim_num} | sma={sma:.1f} inc={np.rad2deg(inc):.1f}deg pattern=({n_planes},{n_planes_sat}) | global={coverage_4fold:.4f} pole={pole_coverage_4fold:.4f}"
    )

    return (i, j, float(coverage_4fold), float(pole_coverage_4fold))


def gridsearch_coverage(
    et0,
    smas,
    incs,
    walker_pattern,
    min_elev_deg=5.0,
    n_user=500,
    dt_sim=60.0,
    pole_lat=-75.0,
    use_hybrid=False,
    parallel=False,
    max_workers=None,
):
    """
    Parallelizes over (sma, inc, walker_pattern) combinations when parallel=True.
    """
    import pylupnt as pnt

    n_sma = len(smas)
    n_incs = len(incs)
    n_res = 2

    sim_num = n_sma * n_incs
    results = np.ones((n_sma, n_incs, n_res), dtype=float) * np.nan

    # Users & ENU matrices (shared, read-only)
    x_user = fibonacci_sphere(n_user) * pnt.R_MOON  # [km]
    lat_users = np.arcsin(x_user[:, 2] / pnt.R_MOON)
    pole_mask = lat_users <= np.deg2rad(pole_lat)
    enu_mats = compute_enu_matrices(x_user)  # (n_users, 3, 3)

    if not parallel:
        # Original serial loop (kept with tqdm)
        for simi in range(sim_num):
            i = simi // (n_incs)
            j = simi % n_incs

            sma = smas[i]
            inc = incs[j]
            pattern = walker_pattern

            _, _, g, p = run_one_coverage_case(
                (
                    i,
                    j,
                    sma,
                    inc,
                    pattern,
                    et0,
                    dt_sim,
                    min_elev_deg,
                    sim_num,
                    simi,
                    x_user,
                    enu_mats,
                    pole_mask,
                    use_hybrid,
                )
            )
            results[i, j, 0] = g
            results[i, j, 1] = p

        return results

    # ---- Parallel path ----
    if max_workers is None:
        max_workers = os.cpu_count() or 2

    # Build task list
    tasks = []
    for simi in range(sim_num):
        i = simi // (n_incs)
        j = simi % n_incs
        tasks.append(
            (
                i,
                j,
                float(smas[i]),
                float(incs[j]),
                tuple(walker_pattern),
                float(et0),
                float(dt_sim),
                float(min_elev_deg),
                sim_num,
                simi,
                x_user,
                enu_mats,
                pole_mask,
                use_hybrid,
            )
        )

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(run_one_coverage_case, t) for t in tasks]
        for fut in as_completed(futures):
            try:
                i, j, g, p = fut.result()
                results[i, j, 0] = g
                results[i, j, 1] = p
            except Exception:
                # On failure, leave NaNs for that case
                pass

    return results


def gridsearch_receiver_noise(et0, smas, incs, P_tx_dBw=15.0, min_elev_deg=5.0, CN0_thresh=30.0):

    pnt.set_lupnt_epoch(0.0)

    # constants
    n_sma = len(smas)
    n_incs = len(incs)
    dt_sim = 60.0

    sim_num = n_sma * n_incs
    results = np.ones((n_sma, n_incs, 6), dtype=float) * np.nan

    # setup user points ----------------------------------------------
    n_user = 250
    x_user = fibonacci_sphere(n_user) * pnt.R_MOON  # [km] user positions on the surface of the Moon
    enu_mats = compute_enu_matrices(x_user)  # (n_users, 3, 3) PA->ENU matrices for user positions

    for i in tqdm(range(n_sma), desc="Gridsearch CN0", unit="sma"):

        sma = smas[i]

        # Skip infeasible eccentricities ------------------------------
        ecc2 = 1 - 5 / 3 * np.cos(incs) ** 2
        ecc = np.sqrt(np.maximum(0, ecc2))
        peri_h = sma * (1 - ecc) - pnt.R_MOON
        valid_idx = np.where((peri_h >= 100.0) * (ecc2 >= 0))[0]
        if len(valid_idx) == 0:
            continue

        n_sat = len(valid_idx)  #

        # deploy satellites ------------------------------------------------
        coes = np.zeros((n_sat, 6))
        coes[:, 0] = sma
        coes[:, 1] = ecc[valid_idx]
        coes[:, 2] = incs[valid_idx]
        coes[:, 3] = 0  # RAAN
        coes[:, 4] = np.deg2rad(-90.0)  # RAAN
        coes[:, 5] = 0.0  # mean anomaly

        T_orbit = 2 * np.pi * np.sqrt(sma**3 / pnt.GM_MOON)
        n_sim = int(T_orbit / dt_sim)
        tspan = np.linspace(0.0, T_orbit, n_sim + 1)
        et = et0 + tspan
        lent = len(tspan)

        # setup Antennas ----------------------------------------------
        lunanet_antennas = [None] * n_sat
        for j in range(n_sat):
            lunanet_antennas[j] = LunanetSatAntenna(P_tx_dBw, coes[j])

        # propagate orbits ------------------------------------------------
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

        rv0_mci = np.zeros((n_sat, 6))
        x_orb = np.zeros((n_sat, lent, 6))
        for si in range(n_sat):
            rv0_op = pnt.classical_to_cart(coes[si, :], pnt.GM_MOON)
            rv0_mci[si] = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)
            x_sat_mci = dynamics.propagate(rv0_mci[si], et)
            x_orb[si] = pnt.convert_frame(et, x_sat_mci, pnt.MOON_CI, pnt.MOON_PA)

        # Compute user-to-satellite vectors
        r_m2u = np.tile(x_user[:, None, None, :], (1, n_sat, lent, 1))  # (n_user, n_sat, lent, 3)
        r_m2s = np.tile(x_orb[None, :, :, :3], (n_user, 1, 1, 1))  # (n_user, n_sat, lent, 3)
        r_u2s = r_m2s - r_m2u  # (n_user, n_sat, lent, 3)

        # Convert to ENU and compute elevation angles ------------------------
        r_u2s_enu = np.matmul(enu_mats[:, None, None, :, :], r_u2s[..., None]).squeeze(-1)
        r_u2s_v = r_u2s_enu.reshape(-1, 3)
        u2s_norm = np.linalg.norm(r_u2s_v, axis=1)
        r_u2s_v /= u2s_norm[:, None]
        r_u2s_enu = r_u2s_v.reshape(n_user, n_sat, lent, 3)

        elev_user = np.arcsin(r_u2s_enu[..., 2])  # elevation = arcsin(z)
        coverage_elev = elev_user >= np.deg2rad(min_elev_deg)  # (n_user, n_sat, lent)

        # compute CN0 ----------------------------------------------------------------
        cn0_links, coverage = compute_lunanet_cn0(
            et,
            x_orb,
            r_u2s,
            elev_user,
            u2s_norm.reshape(n_user, n_sat, lent),
            lunanet_antennas,
            min_elev_deg,
            CN0_thresh,
        )

        # compute ures -----------------------------------------------------------------
        ures_noise = compute_gnss_pseudorange_noise(cn0_links)  # (n_user, n_sat, lent)

        # reshape to (n_sat)

        # for each satellite, compute RMS and 95% ures over all users and time ---------
        ures_noise[~coverage] = np.nan

        # reshape to (n_sat, n_user, lent)
        ures_noise = np.transpose(ures_noise, (1, 0, 2))

        # reshape to (n_sat, n_user*lent)
        ures_noise = ures_noise.reshape(n_sat, -1)

        ure_rms = np.sqrt(np.nanmean(ures_noise**2, axis=1))
        ure_95 = np.nanpercentile(ures_noise, 95, axis=1)
        ure_997 = np.nanpercentile(ures_noise, 99.7, axis=1)

        cn0_links[~coverage_elev] = np.nan
        cn0_links = cn0_links.reshape(n_sat, -1)  # (n_sat, n_user*lent)
        nonnan_num = np.sum(~np.isnan(cn0_links), axis=1)

        cn0_mean = np.nanmean(cn0_links, axis=1)
        cn0_over_20 = np.sum(cn0_links >= 20, axis=1) / nonnan_num
        cn0_over_30 = np.sum(cn0_links >= 30, axis=1) / nonnan_num

        results[i, valid_idx, 0] = ure_rms
        results[i, valid_idx, 1] = ure_95
        results[i, valid_idx, 2] = ure_997
        results[i, valid_idx, 3] = cn0_mean
        results[i, valid_idx, 4] = cn0_over_20
        results[i, valid_idx, 5] = cn0_over_30

    return results


def compute_user_dop(i, lent, nsat_obs, r_u2s_enu, coverage):
    pdop_i = np.full(lent, np.nan)
    gdop_i = np.full(lent, np.nan)

    for j in range(lent):
        if nsat_obs[i, j] < 4:
            continue

        r_enu = r_u2s_enu[i, coverage[i, :, j], j, :]  # (n_visible, 3)
        H = np.zeros((nsat_obs[i, j], 4))
        H[:, 0:3] = -r_enu / np.linalg.norm(r_enu, axis=1)[:, None]
        H[:, 3] = 1.0

        Q = np.linalg.inv(H.T @ H)
        trace3 = np.trace(Q[0:3, 0:3])
        trace4 = np.trace(Q)

        if trace3 > 0:
            pdop_i[j] = np.sqrt(trace3)
        if trace4 > 0:
            gdop_i[j] = np.sqrt(trace4)

    return pdop_i, gdop_i


@njit(parallel=False, fastmath=True, cache=True)
def dop_from_geometry_all_users(r_u2s_enu, coverage, nsat_obs, debug=False):
    """
    r_u2s_enu : (n_user, n_sat, lent, 3) float64, C-contiguous
    coverage  : (n_user, n_sat, lent) bool_
    nsat_obs  : (n_user, lent) int32
    """
    n_user, n_sat, lent, _ = r_u2s_enu.shape
    pdop = np.empty((n_user, lent), dtype=np.float64)
    gdop = np.empty((n_user, lent), dtype=np.float64)
    nan = np.nan

    for i in range(n_user):  # <- single entry/exit, OK to parallelize
        # if debug:
        #     print("computing user {}/{}".format(i+1, n_user))
        for j in range(lent):
            m = nsat_obs[i, j]
            if m < 4:
                pdop[i, j] = nan
                gdop[i, j] = nan
            else:
                # Accumulate H^T H directly to avoid storing H
                HTH = np.zeros((4, 4), dtype=np.float64)
                cnt = 0

                for s in range(n_sat):
                    if coverage[i, s, j]:
                        vx = r_u2s_enu[i, s, j, 0]
                        vy = r_u2s_enu[i, s, j, 1]
                        vz = r_u2s_enu[i, s, j, 2]
                        nrm = math.sqrt(vx * vx + vy * vy + vz * vz)
                        if nrm > 0.0:
                            ux = -vx / nrm
                            uy = -vy / nrm
                            uz = -vz / nrm

                            # outer product r^T r, r = [ux, uy, uz, 1]
                            HTH[0, 0] += ux * ux
                            HTH[0, 1] += ux * uy
                            HTH[0, 2] += ux * uz
                            HTH[0, 3] += ux
                            HTH[1, 0] += uy * ux
                            HTH[1, 1] += uy * uy
                            HTH[1, 2] += uy * uz
                            HTH[1, 3] += uy
                            HTH[2, 0] += uz * ux
                            HTH[2, 1] += uz * uy
                            HTH[2, 2] += uz * uz
                            HTH[2, 3] += uz
                            HTH[3, 0] += ux
                            HTH[3, 1] += uy
                            HTH[3, 2] += uz
                            HTH[3, 3] += 1.0

                            cnt += 1
                # No early break; we just skip extra visibles implicitly.

                if cnt < 4:
                    pdop[i, j] = nan
                    gdop[i, j] = nan
                else:
                    # Tiny Tikhonov regularization avoids singular inverse
                    eps_reg = 1e-10
                    HTH[0, 0] += eps_reg
                    HTH[1, 1] += eps_reg
                    HTH[2, 2] += eps_reg
                    HTH[3, 3] += eps_reg

                    # No try/except here (keeps single-exit structure)
                    Ginv = np.linalg.inv(HTH)
                    t3 = Ginv[0, 0] + Ginv[1, 1] + Ginv[2, 2]
                    t4 = t3 + Ginv[3, 3]
                    pdop[i, j] = math.sqrt(t3) if t3 > 0.0 else nan
                    gdop[i, j] = math.sqrt(t4) if t4 > 0.0 else nan

    return pdop, gdop


def compute_user_dop_numba(r_u2s_enu, coverage, nsat_obs, debug=False):
    """
    Vectorized, jitted replacement for per-user compute_user_dop.
    Returns:
        pdop : (n_user, lent)
        gdop : (n_user, lent)
    """
    # Ensure dtypes/contiguity for Numba
    r_u2s_enu = np.ascontiguousarray(r_u2s_enu, dtype=np.float64)
    coverage = np.ascontiguousarray(coverage, dtype=np.bool_)
    nsat_obs = np.ascontiguousarray(nsat_obs, dtype=np.int32)
    return dop_from_geometry_all_users(r_u2s_enu, coverage, nsat_obs, debug)


def compute_dop_map(
    sma,
    inc,
    walker_pattern,
    et0,
    dt_sim,
    x_user,
    lat_users_vec,
    lon_users_vec,
    enu_mats,
    min_elev_deg=5.0,
    use_hybrid=False,
    parallel=False,
    plot_fig=True,
    sat_fault=False,
    res=None,
    max_lat=None,
    plot_ratio_under=None,
    plot_orbit=True,
    plot_user_view=True,
    debug=False,
):

    pnt.set_lupnt_epoch(0.0)

    # Users & ENU matrices (shared, read-only)
    lat_users = np.arcsin(x_user[:, 2] / pnt.R_MOON)
    lon_users = np.arctan2(x_user[:, 1] / pnt.R_MOON, x_user[:, 0] / pnt.R_MOON)

    if max_lat is not None:
        lat_idx = lat_users <= np.deg2rad(max_lat)
        x_user = x_user[lat_idx]
        enu_mats = enu_mats[lat_idx]

    n_planes, n_planes_sat = walker_pattern

    # Skip infeasible eccentricities
    ecc2 = 1 - 5 / 3 * np.cos(inc) ** 2
    if ecc2 < 0:
        return None, None, None

    ecc = float(np.sqrt(max(0.0, ecc2)))
    peri_h = sma * (1 - ecc) - pnt.R_MOON
    if peri_h < 100.0:
        return None, None, None

    # one orbit period
    T_sid = 27.321661 * 24 * 3600  # [s] one sidereal month
    T_orbit = 2 * np.pi * np.sqrt(sma**3 / pnt.GM_MOON)
    n_orbit = max(1, int(T_sid / T_orbit) + 1)  # simulate at least one full orbit
    n_sim = int(n_orbit * T_orbit / dt_sim)
    tspan = np.linspace(0.0, n_orbit * T_orbit, n_sim + 1)
    et = et0 + tspan
    lent = len(tspan)

    if res is None:

        # deploy satellites
        wsign = 1  # w = 90 deg
        f = 1  # f = 1
        Omega0 = 0.0

        if not use_hybrid:
            x_walker = np.array([sma, ecc, wsign, n_planes, n_planes_sat, f, Omega0], dtype=float)
            coes = setup_walker(x_walker, float_sma=True)

            # if sat_fault:
            #     # remove one satellite from constellation
            #     coes = coes[:-1, :]
        else:
            # combine north and south walker constellations
            wsign = 1
            x_walker_s = np.array([sma, ecc, wsign, n_planes, n_planes_sat, f, Omega0], dtype=float)
            coes_s = setup_walker(x_walker_s, float_sma=True)

            # if sat_fault:
            #     # remove one satellite from constellation
            #     coes_s = coes_s[:-1, :]

            wsign = 0
            x_walker_n = np.array([sma, ecc, wsign, n_planes, n_planes_sat, f, Omega0], dtype=float)
            coes_n = setup_walker(x_walker_n, float_sma=True)

            coes = np.vstack((coes_s, coes_n))

        n_sat = coes.shape[0]

        # dynamics
        dynamics = pnt.NBodyDynamics()
        dynamics.set_integrator(pnt.IntegratorType.RKF45)
        dynamics.set_integrator_params(
            pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10)
        )
        dynamics.add_body(pnt.Body.Moon(2, 2))
        dynamics.add_body(pnt.Body.Earth())
        dynamics.add_body(pnt.Body.Sun())
        dynamics.set_time_step(dt_sim)
        dynamics.set_frame(pnt.MOON_CI)

        # propagate orbits (Moon-CI → Moon-PA)
        if debug:
            print("prpagating orbits...")
        x_orb_pa = np.zeros((n_sat, lent, 6), dtype=float)
        x_orb_mci = np.zeros((n_sat, lent, 6), dtype=float)
        x_orb_pa0 = np.zeros((n_sat, lent, 6), dtype=float)
        for si in range(n_sat):
            rv0_op = pnt.classical_to_cart(coes[si, :], pnt.GM_MOON)
            rv0_mci = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)
            x_orb_mci[si] = dynamics.propagate(rv0_mci, et)
            x_orb_pa[si] = pnt.convert_frame(et, x_orb_mci[si], pnt.MOON_CI, pnt.MOON_PA)
            x_orb_pa0[si] = pnt.convert_frame(
                et0 * np.ones_like(et), x_orb_mci[si], pnt.MOON_CI, pnt.MOON_PA
            )
            if debug:
                print(f" propagated satellite {si+1}/{n_sat}")

        # user-to-satellite vectors
        if debug:
            print("Computing user to satellite vectors...")
        n_user = x_user.shape[0]
        r_m2u = np.tile(x_user[:, None, None, :], (1, n_sat, lent, 1))  # (n_user, n_sat, lent, 3)
        r_m2s = np.tile(x_orb_pa[None, :, :, :3], (n_user, 1, 1, 1))  # (n_user, n_sat, lent, 3)
        r_u2s = r_m2s - r_m2u

        # ENU conversion
        if debug:
            print("computing ENU vectors...")
        r_u2s_enu = np.matmul(enu_mats[:, None, None, :, :], r_u2s[..., None]).squeeze(-1)
        r_u2s_v = r_u2s_enu.reshape(-1, 3)
        u2s_norm = np.linalg.norm(r_u2s_v, axis=1)
        # avoid divide-by-zero if any exact overlaps occur
        u2s_norm[u2s_norm == 0] = 1.0
        r_u2s_v /= u2s_norm[:, None]
        r_u2s_enu = r_u2s_v.reshape(n_user, n_sat, lent, 3)

        elev_user = np.arcsin(r_u2s_enu[..., 2])
        coverage = elev_user >= np.deg2rad(min_elev_deg)
        nsat_obs = np.sum(coverage, axis=1)  # (n_user, lent)

        pdop = np.ones((n_user, lent), dtype=float) * np.nan
        gdop = np.ones((n_user, lent), dtype=float) * np.nan
        if sat_fault:
            pdop = np.ones((n_sat, n_user, lent), dtype=float) * np.nan
            gdop = np.ones((n_sat, n_user, lent), dtype=float) * np.nan
            nsat_obs = np.tile(nsat_obs[None, :, :], (n_sat, 1, 1))  # (n_sat, n_user, lent)

        # compute DOP
        # Parallel execution
        if not parallel:
            # for i in range(n_user):
            #     # pdop[i], gdop[i] = compute_user_dop(i, lent, nsat_obs, r_u2s_enu, coverage)
            if not sat_fault:
                pdop, gdop = compute_user_dop_numba(r_u2s_enu, coverage, nsat_obs, debug)
            else:
                # If simulating a satellite fault, we need to compute DOP with one satellite removed at a time and get the
                for si in range(n_sat):
                    coverage_fault = coverage.copy()
                    coverage_fault[:, si, :] = False
                    nsat_obs[si] = np.sum(coverage_fault, axis=1)
                    pdop[si], gdop[si] = compute_user_dop_numba(
                        r_u2s_enu, coverage_fault, nsat_obs[si], debug
                    )
        else:
            if not sat_fault:
                results = Parallel(n_jobs=-1)(
                    delayed(compute_user_dop)(i, lent, nsat_obs, r_u2s_enu, coverage)
                    for i in tqdm(range(n_user), desc="Computing DOP", unit="user")
                )
                # Collect results
                pdop = np.vstack([res[0] for res in results])
                gdop = np.vstack([res[1] for res in results])
            else:
                # If simulating a satellite fault, we need to compute DOP with one satellite removed at a time and get the
                for si in range(n_sat):
                    coverage_fault = coverage.copy()
                    coverage_fault[:, si, :] = False
                    nsat_obs[si] = np.sum(coverage_fault, axis=1)

                    results = Parallel(n_jobs=-1)(
                        delayed(compute_user_dop)(i, lent, nsat_obs[si], r_u2s_enu, coverage_fault)
                        for i in tqdm(
                            range(n_user), desc=f"Computing DOP with sat {si} fault", unit="user"
                        )
                    )
                    # Collect results
                    for res in results:
                        pdop[si] = np.vstack([res[0] for res in results])
                        gdop[si] = np.vstack([res[1] for res in results])

    else:
        pdop = res["pdop"]
        gdop = res["gdop"]
        nsat_obs = res["nsat_obs"]
        x_orb_mci = res["x_orb_mci"]
        x_orb_pa = res["x_orb_pa"]
        x_orb_pa0 = res["x_orb_pa0"]
        lat_users = res["lat_users"]
        n_sat = x_orb_pa.shape[0]
        lent = x_orb_pa.shape[1]

    if plot_fig and max_lat is None and not sat_fault:
        gdop_p95 = np.nanpercentile(gdop, 95, axis=1)

        # dopmax = min([np.max(pdop_p95[np.isfinite(pdop_p95)]), 100])
        # dopmax = np.ceil(dopmax / 10) * 10
        dopmax = 100
        dopmin = 1.0

        # plot the gdop maps
        x_orb_range = np.linalg.norm(x_orb_pa[:, :, 0:3], axis=2)
        x_orb_lat = np.arcsin(x_orb_pa[:, :, 2] / x_orb_range)
        x_orb_lon = np.arctan2(x_orb_pa[:, :, 1] / x_orb_range, x_orb_pa[:, :, 0] / x_orb_range)

        # Figure 1: GDOP map ----------------------------------------------
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))

        # satellite ground tracks
        dt = tspan[1] - tspan[0]
        orbit_period_idx = int(2 * T_orbit / dt) + 1
        for si in range(n_sat):
            ax.plot(
                np.rad2deg(x_orb_lon[si, :orbit_period_idx:1]),
                np.rad2deg(x_orb_lat[si, :orbit_period_idx:1]),
                "o",
                markersize=0.1,
                color="black",
                alpha=0.8,
            )
        for si in range(n_sat):
            ax.plot(
                np.rad2deg(x_orb_lon[si, 0]),
                np.rad2deg(x_orb_lat[si, 0]),
                "o",
                color="red",
                markersize=3,
            )

        # contour plot
        nlat = len(lat_users_vec)
        nlon = len(lon_users_vec)
        gdop_grid = gdop_p95.reshape((nlon, nlat)).T  # (nlat, nlon)

        if plot_ratio_under is None:  # plot # 95% GDOP values
            p = ax.pcolormesh(
                np.rad2deg(lon_users_vec),
                np.rad2deg(lat_users_vec),
                gdop_grid,
                norm=LogNorm(vmin=dopmin, vmax=dopmax),
                cmap="viridis",
                alpha=0.5,
            )
        else:
            # plot ratio of users under threshold
            ratio_under = np.sum(gdop <= plot_ratio_under, axis=1) / gdop.shape[1]
            ratio_grid = ratio_under.reshape((nlon, nlat)).T
            p = ax.pcolormesh(
                np.rad2deg(lon_users_vec),
                np.rad2deg(lat_users_vec),
                ratio_grid,
                vmin=0.0,
                vmax=1.0,
                cmap="viridis",
                alpha=0.3,
            )
            dopmax = 1.0
            dopmin = 0.0

            # add contour lines for ratio values at 0.25, 0.5, 0.75
            cs = ax.contour(
                np.rad2deg(lon_users_vec),
                np.rad2deg(lat_users_vec),
                ratio_grid,
                levels=[0.5, 0.6, 0.7, 0.8, 0.9],
                # colors="blue",
                cmap="Reds",
                vmin=0.3,
                vmax=1.0,
                linewidths=1.5,
            )
            labels = ax.clabel(cs, inline=True, fontsize=12, fmt="%.2f")
            # for txt in labels:
            #     txt.set_path_effects([
            #         pe.Stroke(linewidth=3, foreground='black'),
            #         pe.Normal()
            #     ])

            # add south pole boundary line (e.g., at -75 deg)
            ax.axhline(y=-75, color="gray", linestyle="--", linewidth=2)

        if use_hybrid:
            ax.set_title(
                "Hybrid Walker pattern: {}x2 plane, {} sat/plane, sma={:.1f}km, inc={:.1f}deg".format(
                    n_planes, n_planes_sat, sma / 1000, np.rad2deg(inc)
                )
            )
        else:
            ax.set_title(
                "Walker pattern: {}x{} sat, sma={:.1f}km, inc={:.1f}deg".format(
                    n_planes, n_planes_sat, sma / 1000, np.rad2deg(inc)
                )
            )
        ax.set_xlabel("Longitude [deg]", fontsize=12)
        ax.set_ylabel("Latitude [deg]", fontsize=12)
        ax.set_xlim([-180, 180])
        ax.set_ylim([-90, 90])
        if plot_ratio_under is None:
            fig.colorbar(p, ax=ax, label="GDOP (95%)")
        else:
            fig.colorbar(p, ax=ax, label=f"Ratio GDOP under {plot_ratio_under}")
        ax.grid(True)
        plt.tight_layout()

        # Figure 2: Satellite Orbits in MCI frame ----------------------------------
        if plot_orbit:
            fig2 = go.Figure()

            # sample only the first two orbit for plotting to avoid overcrowding
            dt = tspan[1] - tspan[0]
            orbit_period_idx = int(T_orbit / dt) + 1
            x_orb_pa0_plot = x_orb_pa0[:, :orbit_period_idx, :]

            if use_hybrid:
                nidx = n_sat // 2
                pnt.plot.plot_orbits(fig2, x_orb_pa0_plot[:nidx], color="blue")
                pnt.plot.plot_orbits(fig2, x_orb_pa0_plot[nidx:], color="orange")
            else:
                if ecc < 0.1:
                    pnt.plot.plot_orbits(fig2, x_orb_pa0_plot, color="green")
                else:
                    pnt.plot.plot_orbits(fig2, x_orb_pa0_plot, color="blue")

            pnt.plot.scatter(fig2, x_orb_pa0_plot[:, 0, :3], color="red")

            pnt.plot.plot_body(
                fig2,
                pnt.MOON,
                size_factor=2,
                alpha=0.5,
            )
            pnt.plot.set_view(fig2, -45, 20, 2.5)
        else:
            fig2 = None

        # Figure 3: Plot satellite view from users at different latitudes ----------------------------------
        lat_plot_deg = [0, -45, -90]  # latitudes to plot
        lon_plot_deg = [0, 0, 0]  # corresponding longitudes to plot (can be adjusted as needed)
        x_user_idxs = []
        x_users_plot = np.zeros((len(lat_plot_deg), 3))
        for u, lat in enumerate(lat_plot_deg):
            x_user_idxs.append(
                np.argmin(
                    np.abs(lat_users - np.deg2rad(lat))
                    + np.abs(lon_users - np.deg2rad(lon_plot_deg[u]))
                )
            )  # find closest user indices to desired latitudes
            x_users_plot[u] = x_user[x_user_idxs[u]]  # (n_users_plot, 3)

            lat_plot_deg[u] = np.rad2deg(lat_users[x_user_idxs[u]])
            lon_plot_deg[u] = np.rad2deg(lon_users[x_user_idxs[u]])

        dt = tspan[1] - tspan[0]
        orbit_period_idx = int(2 * T_orbit / dt) + 1

        if plot_user_view:
            n_users_plot = len(lat_plot_deg)

            # find the x_user indices closest to the desired lat/lon for plotting
            enu_mats_plot = compute_enu_matrices(x_users_plot)  # (n_users_plot, 3, 3)

            fig3 = plt.figure(figsize=(4 * n_users_plot, 8))
            gs = fig3.add_gridspec(2, n_users_plot, height_ratios=[3, 1])
            axes = np.empty((2, n_users_plot), dtype=object)
            for u in range(n_users_plot):
                axes[0, u] = fig3.add_subplot(gs[0, u], projection="polar")
                axes[1, u] = fig3.add_subplot(gs[1, u])

            for u in range(n_users_plot):
                r_m2u_plot = x_users_plot[u]  # (3,)
                r_m2s_plot = x_orb_pa[:, :, :3]  # (n_sat, lent, 3)
                r_u2s_enu_plot = np.zeros_like(r_m2s_plot)  # (n_sat, lent, 3)
                lent = r_m2s_plot.shape[1]
                az_plot = np.zeros((n_sat, lent))
                elev_plot = np.zeros((n_sat, lent))
                for s in range(n_sat):
                    r_u2s_plot = r_m2s_plot[s] - np.tile(r_m2u_plot, (lent, 1))  # (lent, 3)
                    r_u2s_enu_plot[s] = (
                        enu_mats_plot[u] @ r_u2s_plot.T
                    ).T  # (3, lent) = (3, 3) @ (3, lent)
                    los_norms = np.linalg.norm(r_u2s_enu_plot[s], axis=1)  # (lent,)
                    az_plot[s] = np.arctan2(
                        r_u2s_enu_plot[s, :, 0] / los_norms, r_u2s_enu_plot[s, :, 1] / los_norms
                    )  # (lent,)
                    elev_plot[s] = np.arcsin(r_u2s_enu_plot[s, :, 2] / los_norms)  # (lent,)

                # create a polar plot of satellite positions in the sky as seen by the user
                # print(f"User at lat {lat_plot_deg[u]:.1f} deg  long {lon_plot_deg[u]:.1f} deg:")
                # print(" Azimuth (deg):", np.rad2deg(az_plot[:, 0]))
                # print(" Elevation (deg):", np.rad2deg(elev_plot[:, 0]))

                # Plot 1: Skyplot of satellite positions as seen by the user, with azimuth as the angle and elevation as the radius
                ax = axes[0, u]
                # set to polar plot with azimuth as angle and elevation as radius (90-elevation to have 90 deg at center and 0 deg at edge)
                ax.projection = "polar"

                inv_plot = 5  # interval for plotting satellite positions to avoid overcrowding (plot every 5th position)
                len_visible_idxs = len(az_plot[si, :orbit_period_idx:inv_plot])
                visible_idxs = np.zeros((n_sat, len_visible_idxs), dtype=bool)
                colors = [
                    "blue",
                    "orange",
                    "green",
                    "purple",
                    "cyan",
                    "magenta",
                    "yellow",
                    "brown",
                    "pink",
                    "gray",
                    "olive",
                    "teal",
                    "navy",
                    "maroon",
                    "lime",
                    "coral",
                    "gold",
                    "silver",
                    "black",
                    "red",
                ]  # cycle of colors for satellites
                for si in range(n_sat):
                    az_plot_si = az_plot[si, :orbit_period_idx:inv_plot]
                    elev_plot_si = elev_plot[si, :orbit_period_idx:inv_plot]
                    visible_idx = (
                        elev_plot[si, :orbit_period_idx:inv_plot] >= 0
                    )  # only plot visible satellites below
                    ax.plot(
                        az_plot_si[visible_idx],
                        90 - np.rad2deg(elev_plot_si[visible_idx]),
                        "o",
                        markersize=2,
                        alpha=0.2,
                        color=colors[si % len(colors)],
                        label=f"Sat {si}",
                        zorder=1,
                    )
                    visible_idxs[si] = visible_idx

                for si in range(n_sat):
                    if visible_idxs[
                        si, 0
                    ]:  # if the satellite is visible at t=0, plot its position at t=0 in reds
                        ax.scatter(
                            az_plot[si, 0],
                            90 - np.rad2deg(elev_plot[si, 0]),
                            s=50,
                            edgecolors="black",
                            label="t=0",
                            color=colors[si % len(colors)],
                            zorder=10,
                        )

                # plot the initial satellite positions at t=0 in red
                ax.set_title(
                    f"User at lat {lat_plot_deg[u]:.1f} deg  long {lon_plot_deg[u]:.1f} deg",
                    fontsize=12,
                )
                ax.set_ylim(0, 90)
                ax.set_yticks([0, 30, 60, 90])
                ax.set_yticklabels(["90", "60", "30", "0"])
                ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
                ax.set_xticklabels(["0", "45", "90", "135", "180", "225", "270", "315"])
                ax.grid(True)

                # Plot 2: Time series of gdop
                ax = axes[1, u]
                # normal plot of gdop over time for the user
                gdop_user = gdop[x_user_idxs[u], :orbit_period_idx]
                ax.plot(tspan[:orbit_period_idx] / 3600, gdop_user, color="blue")
                ax.set_title(
                    f"GDOP over time for user at lat {lat_plot_deg[u]:.1f} deg", fontsize=12
                )
                ax.set_xlabel("Time [hours]", fontsize=12)
                ax.set_ylabel("GDOP", fontsize=12)
                ax.axhline(y=6, color="red", linestyle="--", label="GDOP=6")
                ax.set_ylim(0, 50)
                ax.grid(True)
                ax.legend()

            plt.tight_layout()
        else:
            fig3 = None

    else:
        fig = None
        fig2 = None
        fig3 = None

    res = {}
    res["pdop"] = pdop
    res["gdop"] = gdop
    res["nsat_obs"] = nsat_obs
    res["x_orb_mci"] = x_orb_mci
    res["x_orb_pa"] = x_orb_pa
    res["x_orb_pa0"] = x_orb_pa0
    res["lat_users"] = lat_users

    return fig, fig2, fig3, res


def gridsearch_dop_polar(sma, incs, n_user, et0, dt_sim, recompute=False, sat_fault=False):

    # precompute user positions and ENU matrices for DOP calculation
    # Users & ENU matrices (shared, read-only)
    lat_users_vec = np.linspace(-np.pi / 2, np.pi / 2, 45)  # dummy, not used
    lon_users_vec = np.linspace(-np.pi, np.pi, 90)  # dummy, not used
    x_user = fibonacci_sphere(n_user) * pnt.R_MOON  # [km] user positions on the surface of the Moon
    #
    enu_mats = compute_enu_matrices(x_user)  # (n_users, 3, 3)

    walker_patterns = []
    # 6 sats
    walker_patterns += [(2, 3), (3, 2)]
    # 8 sats
    walker_patterns += [(2, 4), (4, 2)]
    # 9 sats
    walker_patterns += [(3, 3)]
    # 10 sats
    walker_patterns += [(2, 5), (5, 2)]
    # 12 sats
    walker_patterns += [(3, 4), (4, 3), (2, 6), (6, 2)]
    # 14 sats
    walker_patterns += [(2, 7)]
    # 16 sats
    walker_patterns += [(4, 4), (2, 8), (8, 2)]

    # Compute eccentricity for critical inclination
    use_hybrid = False

    results = {}
    run_sim_pattern = [True] * len(walker_patterns)

    if sat_fault:
        filename = "data/gridsearch_coverage/polar_fault_{0:.0f}.pkl".format(sma * 1e-3)
    else:
        filename = "data/gridsearch_coverage/polar_{0:.0f}.pkl".format(sma * 1e-3)

    if os.path.exists(filename) and not recompute:
        with open(filename, "rb") as f:
            sma, incs, walker_patterns_loaded, results_loaded = pickle.load(f)
        print("Loaded existing results from", filename)

        for wpi, wp in enumerate(walker_patterns):
            if wp not in results_loaded:
                run_sim_pattern[wpi] = True
            else:
                run_sim_pattern[wpi] = False

    for k, walker_pattern in enumerate(walker_patterns):

        print("-------------------------------------")
        print(k, "/", len(walker_patterns), "  pattern:", walker_pattern)
        print("-------------------------------------")

        if not run_sim_pattern[k]:
            results[walker_pattern] = results_loaded[walker_pattern]
            print("Skipping existing pattern:", walker_pattern)
            continue

        results[walker_pattern] = (
            np.zeros((len(incs), 4)) * np.nan
        )  # global coverage, 95% DOP, 95% polar DOP

        for i in range(len(incs)):
            inc = incs[i]
            _, _, _, res = compute_dop_map(
                sma,
                inc,
                walker_pattern,
                et0,
                dt_sim,
                x_user,
                lat_users_vec,
                lon_users_vec,
                enu_mats,
                min_elev_deg=5.0,
                max_lat=-75.0,
                parallel=False,
                use_hybrid=use_hybrid,
                plot_fig=False,
                sat_fault=sat_fault,
            )

            if res is None:
                continue

            # pole coverage
            if sat_fault:
                for si in range(res["pdop"].shape[0]):
                    # coverage
                    pole_cov_si = np.sum(res["nsat_obs"][si] >= 4) / res["nsat_obs"][si].size
                    if si == 0:
                        pole_cov = pole_cov_si
                    else:
                        pole_cov = min(pole_cov, pole_cov_si)

                    # 95 percentile DOP
                    pole_dop_p95_si = np.nanpercentile(
                        res["pdop"][si].flatten(), 95
                    )  # 95 percentile DOP for each user
                    if si == 0:
                        pole_dop_p95 = pole_dop_p95_si
                    else:
                        pole_dop_p95 = max(pole_dop_p95, pole_dop_p95_si)

                    # GDOP <= 6 ratio
                    pole_dop_under6_si = np.sum(res["gdop"][si] <= 6) / res["gdop"][si].size
                    if si == 0:
                        pole_dop_under6 = pole_dop_under6_si
                    else:
                        pole_dop_under6 = min(pole_dop_under6, pole_dop_under6_si)

                    # GDOP <= 3 ratio
                    pole_dop_under3_si = np.sum(res["gdop"][si] <= 3) / res["gdop"][si].size
                    if si == 0:
                        pole_dop_under3 = pole_dop_under3_si
                    else:
                        pole_dop_under3 = min(pole_dop_under3, pole_dop_under3_si)

            else:
                pole_cov = np.sum(res["nsat_obs"] >= 4) / res["nsat_obs"].size

                # 95 percentile DOP
                pole_dop_p95 = np.nanpercentile(
                    res["pdop"].flatten(), 95
                )  # 95 percentile DOP for each user

                # GDOP <= 6 ratio
                pole_dop_under6 = np.sum(res["gdop"] <= 6) / res["gdop"].size

                # GDOP <= 3 ratio
                pole_dop_under3 = np.sum(res["gdop"] <= 3) / res["gdop"].size

            results[walker_pattern][i, 0] = pole_cov
            results[walker_pattern][i, 1] = pole_dop_p95
            results[walker_pattern][i, 2] = pole_dop_under6
            results[walker_pattern][i, 3] = pole_dop_under3

            print(
                "sma: {:.1f} inc: {:.1f} | pole_cov: {:.3f} pole_dop_p95: {:.1f} pole_dop_under6: {:.3f} | pole_dop_under3: {:.3f}".format(
                    sma,
                    np.rad2deg(inc),
                    pole_cov,
                    pole_dop_p95,
                    pole_dop_under6,
                    pole_dop_under3,
                )
            )

    # Save results
    with open(filename, "wb") as f:
        pickle.dump((sma, incs, walker_patterns, results), f)

    return sma, incs, walker_patterns, results


def gridsearch_dop_hybrid(sma, incs, n_user, et0, dt_sim, recompute=False, sat_fault=False):
    # precompute user positions and ENU matrices for DOP calculation
    # Users & ENU matrices (shared, read-only)
    lat_users_vec = np.linspace(-np.pi / 2, np.pi / 2, 45)  # dummy, not used
    lon_users_vec = np.linspace(-np.pi, np.pi, 90)  # dummy, not used
    x_user = fibonacci_sphere(n_user) * pnt.R_MOON  # [km] user positions on the surface of the Moon
    #
    enu_mats = compute_enu_matrices(x_user)  # (n_users, 3, 3)

    # hybrid Walker patterns to test
    walker_patterns = [(2, 4), (4, 2), (3, 3), (2, 5), (5, 2), (3, 4), (4, 3), (2, 6), (6, 2)]

    # Compute eccentricity for critical inclination
    use_hybrid = True

    results = {}
    run_sim_pattern = [True] * len(walker_patterns)

    if sat_fault:
        filename = "data/gridsearch_coverage/hybrid_fault_{0:.0f}.pkl".format(sma * 1e-3)
    else:
        filename = "data/gridsearch_coverage/hybrid_{0:.0f}.pkl".format(sma * 1e-3)

    if os.path.exists(filename) and not recompute:
        with open(filename, "rb") as f:
            sma, incs, walker_patterns_loaded, results_loaded = pickle.load(f)
        print("Loaded existing results from", filename)

        for wpi, wp in enumerate(walker_patterns):
            if wp not in results_loaded:
                run_sim_pattern[wpi] = True
            else:
                run_sim_pattern[wpi] = False

    for k, walker_pattern in enumerate(walker_patterns):

        print("-------------------------------------")
        print(k, "/", len(walker_patterns), "  pattern:", walker_pattern)
        print("-------------------------------------")

        if not run_sim_pattern[k]:
            results[walker_pattern] = results_loaded[walker_pattern]
            print("Skipping existing pattern:", walker_pattern)
            continue

        results[walker_pattern] = (
            np.zeros((len(incs), 8)) * np.nan
        )  # global coverage, 95% DOP, 95% polar DOP

        for i, inc in enumerate(incs):
            print("Incination: {:.1f} deg".format(np.rad2deg(inc)))
            _, _, _, res = compute_dop_map(
                sma,
                inc,
                walker_pattern,
                et0,
                dt_sim,
                x_user,
                lat_users_vec,
                lon_users_vec,
                enu_mats,
                min_elev_deg=5.0,
                use_hybrid=use_hybrid,
                plot_fig=False,
                sat_fault=sat_fault,
                parallel=False,
            )

            if res is None:
                continue

            if sat_fault:
                lent = res["gdop"][0].shape[1]  # number of users for coverage calculations
            else:
                lent = res["gdop"].shape[1]

            # 4-fold global coverage
            if sat_fault:
                # If simulating a satellite fault, we need to compute coverage with one satellite removed at a time and get the worst-case coverage
                for si in range(res["pdop"].shape[0]):
                    # coverage
                    glob_cov_si = np.sum(res["nsat_obs"][si] >= 4) / res["nsat_obs"][si].size
                    if si == 0:
                        glob_cov = glob_cov_si
                    else:
                        glob_cov = min(glob_cov, glob_cov_si)

                # 95 percentile DOP                glob_dop_95 = 0.0
                glob_dop_95 = 0.0
                for si in range(res["pdop"].shape[0]):
                    glob_dop_95_si = np.nanpercentile(
                        res["pdop"][si].flatten(), 95
                    )  # 95 percentile DOP for each user
                    if si == 0:
                        glob_dop_95 = glob_dop_95_si
                    else:
                        glob_dop_95 = max(glob_dop_95, glob_dop_95_si)

                # GDOP <= 6 ratio
                glob_dop_under6 = 1.0
                for si in range(res["gdop"].shape[0]):
                    glob_dop_under6_si = np.sum(res["gdop"][si] <= 6) / res["gdop"][si].size
                    if si == 0:
                        glob_dop_under6 = glob_dop_under6_si
                    else:
                        glob_dop_under6 = min(glob_dop_under6, glob_dop_under6_si)

                # GDOP <= 3 ratio
                glob_dop_under3 = 1.0
                for si in range(res["gdop"].shape[0]):
                    glob_dop_under3_si = np.sum(res["gdop"][si] <= 3) / res["gdop"][si].size
                    if si == 0:
                        glob_dop_under3 = glob_dop_under3_si
                    else:
                        glob_dop_under3 = min(glob_dop_under3, glob_dop_under3_si)

                # 95% polar DOP
                lat_idx = res["lat_users"] <= np.deg2rad(-75)  # users near the south pole
                pole_dop_95 = 0.0
                for si in range(res["gdop"].shape[0]):
                    pole_dop_95_si = np.nanpercentile(
                        res["gdop"][si][lat_idx].flatten(), 95
                    )  # 95 percentile DOP for each user
                    if si == 0:
                        pole_dop_95 = pole_dop_95_si
                    else:
                        pole_dop_95 = max(pole_dop_95, pole_dop_95_si)

                # pole coverage
                pole_cov = 1.0
                for si in range(res["nsat_obs"].shape[0]):
                    pole_cov_si = np.sum(res["nsat_obs"][si][lat_idx] >= 4) / np.sum(lat_idx) / lent
                    if si == 0:
                        pole_cov = pole_cov_si
                    else:
                        pole_cov = min(pole_cov, pole_cov_si)

                # pole gdop <= 6 ratio
                pole_dop_under6 = 1.0
                for si in range(res["gdop"].shape[0]):
                    pole_dop_under6_si = (
                        np.sum(res["gdop"][si][lat_idx] <= 6) / np.sum(lat_idx) / lent
                    )
                    if si == 0:
                        pole_dop_under6 = pole_dop_under6_si
                    else:
                        pole_dop_under6 = min(pole_dop_under6, pole_dop_under6_si)

                # pole gdop <= 3 ratio
                pole_dop_under3 = 1.0
                for si in range(res["gdop"].shape[0]):
                    pole_dop_under3_si = (
                        np.sum(res["gdop"][si][lat_idx] <= 3) / np.sum(lat_idx) / lent
                    )
                    if si == 0:
                        pole_dop_under3 = pole_dop_under3_si
                    else:
                        pole_dop_under3 = min(pole_dop_under3, pole_dop_under3_si)

            else:
                glob_cov = np.sum(res["nsat_obs"] >= 4) / res["nsat_obs"].size

                # 95 percentile DOP
                glob_dop_95 = np.nanpercentile(
                    res["pdop"].flatten(), 95
                )  # 95 percentile DOP for each user

                # GDOP <= 6 ratio
                glob_dop_under6 = np.sum(res["gdop"] <= 6) / res["gdop"].size

                # GDOP <= 3 ratio
                glob_dop_under3 = np.sum(res["gdop"] <= 3) / res["gdop"].size

                # 95% polar DOP
                lat_idx = res["lat_users"] <= np.deg2rad(-75)  # users near the south pole
                pole_dop_95 = np.nanpercentile(
                    res["gdop"][lat_idx].flatten(), 95
                )  # 95 percentile DOP for each user

                # pole coverage
                pole_cov = np.sum(res["nsat_obs"][lat_idx] >= 4) / np.sum(lat_idx) / lent

                # pole gdop <= 6 ratio
                pole_dop_under6 = np.sum(res["gdop"][lat_idx] <= 6) / np.sum(lat_idx) / lent

                # pole gdop <= 3 ratio
                pole_dop_under3 = np.sum(res["gdop"][lat_idx] <= 3) / np.sum(lat_idx) / lent

            results[walker_pattern][i, 0] = glob_cov
            results[walker_pattern][i, 1] = glob_dop_95
            results[walker_pattern][i, 2] = glob_dop_under6
            results[walker_pattern][i, 3] = glob_dop_under3
            results[walker_pattern][i, 4] = pole_dop_95
            results[walker_pattern][i, 5] = pole_cov
            results[walker_pattern][i, 6] = pole_dop_under6
            results[walker_pattern][i, 7] = pole_dop_under3

            print(
                "sma: {:.1f} inc: {:.1f} | glob_cov: {:.3f} glob_dop95: {:.1f} glob_dop<=6: {:.3f} | pole_dop95: {:.1f} pole_cov: {:.3f} pole_dop<=6: {:.3f}".format(
                    sma,
                    np.rad2deg(inc),
                    glob_cov,
                    glob_dop_95,
                    glob_dop_under6,
                    pole_dop_95,
                    pole_cov,
                    pole_dop_under6,
                )
            )

    # Save results
    with open(filename, "wb") as f:
        pickle.dump((sma, incs, walker_patterns, results), f)

    return sma, incs, walker_patterns, results


def gridsearch_dop_circular(smas, n_user, et0, dt_sim, recompute=False, sat_fault=False):

    # precompute user positions and ENU matrices for DOP calculation
    # Users & ENU matrices (shared, read-only)
    lat_users_vec = np.linspace(-np.pi / 2, np.pi / 2, 45)  # dummy, not used
    lon_users_vec = np.linspace(-np.pi, np.pi, 90)  # dummy, not used
    x_user = fibonacci_sphere(n_user) * pnt.R_MOON  # [km] user positions on the surface of the Moon
    #
    enu_mats = compute_enu_matrices(x_user)  # (n_users, 3, 3)

    walker_patterns = []
    # 8 sats
    walker_patterns += [(2, 4), (4, 2)]
    # 9 sats
    walker_patterns += [(3, 3)]
    # 10 sats
    walker_patterns += [(2, 5), (5, 2)]
    # 12 sats
    walker_patterns += [(3, 4), (4, 3), (2, 6), (6, 2)]
    # 16 sats
    walker_patterns += [(4, 4)]
    # 18 sats
    walker_patterns += [(6, 3), (3, 6)]
    # 20 sats
    walker_patterns += [(4, 5), (5, 4)]
    # 21 sats
    walker_patterns += [(7, 3), (3, 7)]
    # 24 sats
    walker_patterns += [(8, 3), (3, 8), (4, 6), (6, 4)]
    # 25 sats
    walker_patterns += [(5, 5)]

    # Compute eccentricity for critical inclination
    eps = 1e-6
    inc = np.arccos(np.sqrt(3 / 5)) + eps  # critical inclination
    ecc = np.sqrt(1 - 5 / 3 * np.cos(inc) ** 2)
    print("Eccentricity for inc={0:.1f} deg: {1:.4f}".format(np.rad2deg(inc), ecc))
    use_hybrid = False

    results = {}
    run_sim_pattern = [True] * len(walker_patterns)

    if sat_fault:
        filename = "data/gridsearch_coverage/circular_fault.pkl"
    else:
        filename = "data/gridsearch_coverage/circular.pkl"

    if os.path.exists(filename) and not recompute:
        with open(filename, "rb") as f:
            smas, inc, walker_patterns_loaded, results_loaded = pickle.load(f)
        print("Loaded existing results from", filename)

        for wpi, wp in enumerate(walker_patterns):
            if wp not in results_loaded:
                run_sim_pattern[wpi] = True
            else:
                run_sim_pattern[wpi] = False

    for k, walker_pattern in enumerate(walker_patterns):

        print("-------------------------------------")
        print(k, "/", len(walker_patterns), "  pattern:", walker_pattern)
        print("-------------------------------------")

        if not run_sim_pattern[k]:
            results[walker_pattern] = results_loaded[walker_pattern]
            print("Skipping existing pattern:", walker_pattern)
            continue

        results[walker_pattern] = (
            np.zeros((len(smas), 8)) * np.nan
        )  # global coverage, 95% DOP, 95% polar DOP

        for i, sma in enumerate(smas):
            _, _, _, res = compute_dop_map(
                sma,
                inc,
                walker_pattern,
                et0,
                dt_sim,
                x_user,
                lat_users_vec,
                lon_users_vec,
                enu_mats,
                min_elev_deg=5.0,
                use_hybrid=use_hybrid,
                plot_fig=False,
                parallel=False,
                sat_fault=sat_fault,
            )

            if res is None:
                continue

            if sat_fault:
                lent = res["gdop"][0].shape[1]  # number of users for coverage calculations
            else:
                lent = res["gdop"].shape[1]

            if sat_fault:
                # If simulating a satellite fault, we need to compute coverage with one satellite removed at a time and get the worst-case coverage
                for si in range(res["gdop"].shape[0]):
                    # coverage
                    glob_cov_si = np.sum(res["nsat_obs"][si] >= 4) / res["nsat_obs"][si].size
                    if si == 0:
                        glob_cov = glob_cov_si
                    else:
                        glob_cov = min(glob_cov, glob_cov_si)

                # 95 percentile DOP                glob_dop_95 = 0.0
                glob_dop_95 = 0.0
                for si in range(res["gdop"].shape[0]):
                    glob_dop_95_si = np.nanpercentile(
                        res["gdop"][si].flatten(), 95
                    )  # 95 percentile DOP for each user
                    if si == 0:
                        glob_dop_95 = glob_dop_95_si
                    else:
                        glob_dop_95 = max(glob_dop_95, glob_dop_95_si)

                # GDOP <= 6 ratio
                glob_dop_under6 = 1.0
                for si in range(res["gdop"].shape[0]):
                    glob_dop_under6_si = np.sum(res["gdop"][si] <= 6) / res["gdop"][si].size
                    if si == 0:
                        glob_dop_under6 = glob_dop_under6_si
                    else:
                        glob_dop_under6 = min(glob_dop_under6, glob_dop_under6_si)

                # GDOP <= 3 ratio
                glob_dop_under3 = 1.0
                for si in range(res["gdop"].shape[0]):
                    glob_dop_under3_si = np.sum(res["gdop"][si] <= 3) / res["gdop"][si].size
                    if si == 0:
                        glob_dop_under3 = glob_dop_under3_si
                    else:
                        glob_dop_under3 = min(glob_dop_under3, glob_dop_under3_si)

                # 95% polar DOP
                lat_idx = res["lat_users"] <= np.deg2rad(-75)  # users near the south pole
                pole_dop_95 = 0.0
                for si in range(res["gdop"].shape[0]):
                    pole_dop_95_si = np.nanpercentile(
                        res["gdop"][si][lat_idx].flatten(), 95
                    )  # 95 percentile DOP for each user
                    if si == 0:
                        pole_dop_95 = pole_dop_95_si
                    else:
                        pole_dop_95 = max(pole_dop_95, pole_dop_95_si)

                # pole coverage
                pole_cov = 1.0
                for si in range(res["nsat_obs"].shape[0]):
                    pole_cov_si = np.sum(res["nsat_obs"][si][lat_idx] >= 4) / np.sum(lat_idx) / lent
                    if si == 0:
                        pole_cov = pole_cov_si
                    else:
                        pole_cov = min(pole_cov, pole_cov_si)

                # pole gdop <= 6 ratio
                pole_dop_under6 = 1.0
                for si in range(res["gdop"].shape[0]):
                    pole_dop_under6_si = (
                        np.sum(res["gdop"][si][lat_idx] <= 6) / np.sum(lat_idx) / lent
                    )
                    if si == 0:
                        pole_dop_under6 = pole_dop_under6_si
                    else:
                        pole_dop_under6 = min(pole_dop_under6, pole_dop_under6_si)

                # pole gdop <= 3 ratio
                pole_dop_under3 = 1.0
                for si in range(res["gdop"].shape[0]):
                    pole_dop_under3_si = (
                        np.sum(res["gdop"][si][lat_idx] <= 3) / np.sum(lat_idx) / lent
                    )
                    if si == 0:
                        pole_dop_under3 = pole_dop_under3_si
                    else:
                        pole_dop_under3 = min(pole_dop_under3, pole_dop_under3_si)

            else:
                # 4-fold global coverage
                glob_cov = np.sum(res["nsat_obs"] >= 4) / res["nsat_obs"].size

                # 95 percentile DOP
                glob_dop_95 = np.nanpercentile(
                    res["gdop"].flatten(), 95
                )  # 95 percentile DOP for each user

                # GDOP <= 6 ratio
                glob_dop_under6 = np.sum(res["gdop"] <= 6) / res["gdop"].size

                # GDOP <= 3 ratio
                glob_dop_under3 = np.sum(res["gdop"] <= 3) / res["gdop"].size

                # 95% polar DOP
                lat_idx = res["lat_users"] <= np.deg2rad(-75)  # users near the south pole
                pole_dop_95 = np.nanpercentile(
                    res["gdop"][lat_idx].flatten(), 95
                )  # 95 percentile DOP for each user

                # pole coverage
                pole_cov = np.sum(res["nsat_obs"][lat_idx] >= 4) / np.sum(lat_idx) / lent

                # pole gdop <= 6 ratio
                pole_dop_under6 = np.sum(res["gdop"][lat_idx] <= 6) / np.sum(lat_idx) / lent

                # pole gdop <= 3 ratio
                pole_dop_under3 = np.sum(res["gdop"][lat_idx] <= 3) / np.sum(lat_idx) / lent

            results[walker_pattern][i, 0] = glob_cov
            results[walker_pattern][i, 1] = glob_dop_95
            results[walker_pattern][i, 2] = glob_dop_under6
            results[walker_pattern][i, 3] = glob_dop_under3
            results[walker_pattern][i, 4] = pole_dop_95
            results[walker_pattern][i, 5] = pole_cov
            results[walker_pattern][i, 6] = pole_dop_under6
            results[walker_pattern][i, 7] = pole_dop_under3

            print(
                "sma: {:.1f} km inc: {:.1f} | glob_cov: {:.3f} glob_dop95: {:.1f} glob_dop<=6: {:.3f} | pole_dop95: {:.1f} pole_cov: {:.3f} pole_dop<=6: {:.3f}".format(
                    sma / 1000,
                    np.rad2deg(inc),
                    glob_cov,
                    glob_dop_95,
                    glob_dop_under6,
                    pole_dop_95,
                    pole_cov,
                    pole_dop_under6,
                )
            )

    # Save results
    with open(filename, "wb") as f:
        pickle.dump((smas, inc, walker_patterns, results), f)

    return smas, inc, walker_patterns, results
