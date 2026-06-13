import pylupnt as pnt
import numpy as np
import matplotlib.pyplot as plt


def eci_to_serot(t_epoch, x_eci):
    """
    Convert ECI state to Sun-Earth fixed frame.

    Parameters
    ----------
    t_epoch : datetime [1, N]
        Epoch time in UTC.
    x_eci : ndarray [N, 6]
        ECI state vector [x, y, z, vx, vy, vz] in km and km/s.
    Returns
    -------
    x_sef : ndarray [N, 6]
        State vector in Sun-Earth fixed frame [x, y, z, vx, vy, vz] in km and km/s.
    """
    sun_rv = pnt.get_body_pos_vel(t_epoch, pnt.EARTH, pnt.SUN, pnt.ECI)  # (N, 6)

    assert x_eci.shape[1] == 6 and sun_rv.shape[1] == 6
    assert x_eci.shape[0] == sun_rv.shape[0]
    N = x_eci.shape[0]

    r_sc = x_eci[:, 0:3]
    v_sc = x_eci[:, 3:6]

    r_sun = sun_rv[:, 0:3]
    v_sun = sun_rv[:, 3:6]

    # Unit vectors for the SEF frame
    r_sun_norm = np.linalg.norm(r_sun, axis=1)
    if np.any(r_sun_norm < 1e-12):
        raise ValueError("Sun position has near-zero norm for some samples.")

    x_hat = -r_sun / r_sun_norm[:, None]  # +x points away from the Sun

    h_sun = np.cross(r_sun, v_sun)
    h_norm = np.linalg.norm(h_sun, axis=1)
    if np.any(h_norm < 1e-15):
        raise ValueError("Sun angular momentum near zero; cannot define z-axis.")

    z_hat = h_sun / h_norm[:, None]
    y_hat = np.cross(z_hat, x_hat)

    # Re-orthonormalize y_hat (and fix tiny numerical drift)
    y_norm = np.linalg.norm(y_hat, axis=1)
    if np.any(y_norm < 1e-15):
        raise ValueError("Degenerate basis: x_hat parallel to z_hat for some samples.")
    y_hat /= y_norm[:, None]

    # Recompute z_hat to ensure orthonormal triad
    z_hat = np.cross(x_hat, y_hat)

    # Rotation matrix R such that v_SEF = R * v_ECI
    # (stack basis vectors as rows)
    R = np.stack([x_hat, y_hat, z_hat], axis=1)  # shape (N, 3, 3)

    # Angular velocity of the rotating SEF frame (in ECI components)
    omega = h_sun / (r_sun_norm**2)[:, None]  # shape (N, 3)

    # Transport term for velocities: v_rot = v_eci - omega × r_eci
    omega_cross_r = np.cross(omega, r_sc)
    v_rel = v_sc - omega_cross_r

    # Apply rotation
    r_sef = np.einsum("nij,nj->ni", R, r_sc)
    v_sef = np.einsum("nij,nj->ni", R, v_rel)

    out = np.empty_like(x_eci)
    out[:, 0:3] = r_sef
    out[:, 3:6] = v_sef
    return out


def propagate_patched_conics(et, rv0_mci, debug=False):
    """
    Propagate a low-energy transfer trajectory from Earth to Moon using patched conics.
    Once it enters the Moon's sphere of influence, it switches to a two-body propagation around the Moon.
    Otherwise, it propagates around the Earth considering the Sun and Earth perturbation.
    """
    et0 = et[0]
    etf = et[-1]
    lent = len(et)
    et_query = np.linspace(
        etf - 100, et0 + 100, 100 * lent
    )  # ephemeris times for Sun and Earth positions
    earth_rv = pnt.get_body_pos_vel(et, pnt.MOON, pnt.EARTH, pnt.MOON_CI)
    moon_rv = pnt.get_body_pos_vel(et, pnt.EARTH, pnt.MOON, pnt.ECI)
    r_em = np.linalg.norm(earth_rv[:, 0:3], axis=1)
    a_MOON = np.mean(r_em)  # km  # semi-major axis of the Moon orbit
    r_SOI_MOON = a_MOON * (pnt.GM_MOON / pnt.GM_EARTH) ** (2 / 5)  # km

    def hit_or_too_far(t, r_moon, r_earth):
        # cond1: pass too close to Earth
        if r_earth < (pnt.R_EARTH + 200e3):  # 200 km altitude from Earth center
            if debug:
                print(f" Hit the Earth: r_earth={r_earth/1000:.1f} km")
            return True  # terminate when getting too close to Earth

        # cond2: pass too close to Moon after 5 days
        if (r_moon < (pnt.R_MOON + 500e3)) and (
            (t - et[0]) > 5 * 24 * 3600
        ):  # 500 km altitude from Moon center
            if debug:
                print(f" Hit the Moon: r_moon={r_moon/1000:.1f} km after 5 days")
            return True

        # cond3: hit the Moon
        if r_moon < (pnt.R_MOON + 10e3):
            if debug:
                print(f" Hit the Moon: r_moon={r_moon/1000:.1f} km")
            return True  # terminate when hitting the Moon

        # cond4: go too far from both Earth and Moon
        if r_earth > (a_MOON * 5):
            if debug:
                print(f" Went too far from Earth: r_earth={r_earth/1000:.1f} km")
            return True  # terminate when going too far from both Earth and Moon

        return False

    # Termination condition function
    def terminate_func(t, x, center="Moon", debug=debug):
        if center == "Moon":
            r_moon = np.linalg.norm(x[0:3])
            earth_rv_mci = pnt.get_body_pos_vel(t, pnt.MOON, pnt.EARTH, pnt.MOON_CI)
            r_earth = np.linalg.norm(
                x[:3] - earth_rv_mci[:3]
            )  # relative to Earth (s/c - moon) - (Earth - Moon)
        elif center == "Earth":
            r_earth = np.linalg.norm(x[0:3])
            moon_rv_eci = pnt.get_body_pos_vel(t, pnt.EARTH, pnt.MOON, pnt.ECI)
            r_moon = np.linalg.norm(
                x[:3] - moon_rv_eci[:3]
            )  # relative to Moon (s/c - Earth) - (Moon - Earth)
        else:
            raise ValueError("center must be 'Moon' or 'Earth'")

        # switch patched conics if entering Moon's SOI
        if r_moon <= (r_SOI_MOON - 1000e3) and center == "Earth":
            if debug:
                print(f" Entered Moon's SOI at t={t}, r_moon={r_moon:.1f} km")
            return True  # terminate to switch to Moon-centered propagation

        if r_moon > (r_SOI_MOON + 1000e3) and center == "Moon":
            if debug:
                print(f" Exited Moon's SOI at t={t}, r_moon={r_moon:.1f} km")
            return True  # terminate to switch to Earth-centered propagation

        if hit_or_too_far(t, r_moon, r_earth):
            return True  # terminate if hitting the Earth or Moon surface

        return False  # terminate when reaching 1880 km from Moon center

    iparam_e = pnt.IntegratorParams(max_iter=1000, abstol=1e-10, reltol=1e-10)
    iparam_m = pnt.IntegratorParams(max_iter=1000, abstol=1e-10, reltol=1e-10)
    iparam_e.set_terminate_if(
        lambda t, x: terminate_func(t, x, center="Earth", debug=debug)
    )  # terminate if altitude < 200 km from Earth or < 500 km from Moon after 5 days
    iparam_m.set_terminate_if(
        lambda t, x: terminate_func(t, x, center="Moon", debug=debug)
    )  # terminate if altitude < 500 km from Moon

    # dynamics
    dynamics_e = pnt.NBodyDynamics()
    dynamics_e.set_integrator(pnt.IntegratorType.RKF45)
    dynamics_e.set_integrator_params(iparam_e)
    dynamics_e.add_body(pnt.Body.Earth())
    dynamics_e.add_body(pnt.Body.Moon())
    dynamics_e.add_body(pnt.Body.Sun())
    dynamics_e.set_time_step(300)
    dynamics_e.set_frame(pnt.ECI)

    dynamics_m = pnt.NBodyDynamics()
    dynamics_m.set_integrator(pnt.IntegratorType.RKF45)
    dynamics_m.set_integrator_params(iparam_m)
    dynamics_m.add_body(pnt.Body.Moon())
    dynamics_m.add_body(pnt.Body.Earth())
    dynamics_m.add_body(pnt.Body.Sun())
    dynamics_m.set_time_step(300)
    dynamics_m.set_frame(pnt.MOON_CI)

    # propagate backward from the lunar orbit insertion point
    terminate = False
    terminate_fail = False
    curr_dyn = dynamics_m  # start with Moon-centered propagation
    rv0 = rv0_mci
    curr_center = "Moon"

    x_eci = np.array([])  # propagated states in ECI
    x_mci = np.array([])  # propagated states in MCI

    et0_prop = et0
    et_prop = et.copy()

    while not terminate:
        if curr_center == "Moon":
            xprop_bw_mci, info = curr_dyn.propagate_with_info(rv0, et0_prop, et_prop)
            n_prop = xprop_bw_mci.shape[0]
            earth_rv_mci = pnt.get_body_pos_vel(et_prop[:n_prop], pnt.MOON, pnt.EARTH, pnt.MOON_CI)
            xprop_bw_eci = (
                xprop_bw_mci - earth_rv_mci
            )  # relative to Earth (x - MOON) - (Earth - MOON)
        else:
            xprop_bw_eci, info = curr_dyn.propagate_with_info(rv0, et0_prop, et_prop)
            n_prop = xprop_bw_eci.shape[0]
            moon_rv_eci = pnt.get_body_pos_vel(et_prop[:n_prop], pnt.EARTH, pnt.MOON, pnt.ECI)
            xprop_bw_mci = (
                xprop_bw_eci - moon_rv_eci
            )  # relative to Moon (x - Earth) - (Moon - Earth)

        r_moon_final = np.linalg.norm(xprop_bw_mci[-1, 0:3])
        r_earth_final = np.linalg.norm(xprop_bw_eci[-1, 0:3])

        if debug:
            print(xprop_bw_mci[-1, :])
            print(
                f"  Propagated {n_prop} steps | center: {curr_center} | t: {et_prop[n_prop-1]} | r_moon: {r_moon_final:.1f} km | r_earth: {r_earth_final:.1f} km"
            )

        if info.reason == pnt.TerminationReason.ReachedTf:
            terminate = True  # finished propagation
        else:
            # if the termination is due to patched conics switch, update the dynamics and initial state
            if not hit_or_too_far(et_prop[n_prop - 1], r_moon_final, r_earth_final):
                if curr_center == "Moon":
                    curr_center = "Earth"
                    curr_dyn = dynamics_e
                    rv0 = xprop_bw_eci[-1, :]  # relative to Earth
                else:
                    curr_center = "Moon"
                    curr_dyn = dynamics_m
                    rv0 = xprop_bw_mci[-1, :]  # relative to Moon
                # update the epoch as well
                et0_prop = et_prop[n_prop - 1]
                et_prop = et_prop[n_prop - 1 :]
                if debug:
                    print("   et_prop length: ", len(et_prop))
                terminate = False  # continue propagation
            else:
                terminate = (
                    True  # finished propagation due to other reasons (e.g., hit the Moon or Earth)
                )
                terminate_fail = True  # mark as failed propagation
                if debug:
                    print(
                        f"  Terminate propagation at t={et_prop[n_prop-1]} | r_moon={r_moon_final:.1f} km due to {info.reason}"
                    )

        # exntend the results
        x_eci = np.vstack((x_eci, xprop_bw_eci[1:])) if x_eci.size else xprop_bw_eci
        x_mci = np.vstack((x_mci, xprop_bw_mci[1:])) if x_mci.size else xprop_bw_mci

    return x_eci, x_mci, terminate_fail


from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from itertools import repeat
import numpy as np
import csv
import os


# import pandas as pd  # optional if you want Parquet
def k_step_minima_idxs(r, k=2, eps=0.0):
    """
    Return indices c in r where r strictly decreases for k steps BEFORE c
    and strictly increases for k steps AFTER c:

        r[c-k] > r[c-k+1] > ... > r[c-1] < r[c] < r[c+1] < ... < r[c+k]

    Parameters
    ----------
    r : array_like
        1D sequence of values (e.g., radius vs time).
    k : int
        Number of consecutive steps on each side to enforce (k >= 1).
    eps : float
        Strictness tolerance. Use small positive eps to ignore tiny noise:
        decreases require diff < -eps; increases require diff > eps.

    Returns
    -------
    idx : np.ndarray (int)
        Indices into r of all k-step local minima.
    """
    r = np.asarray(r)
    n = r.size
    if n < 2 * k + 1:
        # compute for k that fits
        k = (n - 1) // 2
        return k_step_minima_idxs(r, k=k, eps=eps)

    if k < 1 or n < 2 * k + 1:
        return np.array([], dtype=int)

    d = np.diff(r)  # length n-1

    # Build "has k consecutive negatives/positives" without sliding_window_view
    L = (n - 1) - k + 1  # = n - k
    decr_k = np.ones(L, dtype=bool)
    incr_k = np.ones(L, dtype=bool)
    for j in range(k):
        decr_k &= d[j : j + L] < -eps  # windows starting at j..j+L-1
        incr_k &= d[j : j + L] > eps

    # Candidate centers c run from k .. (n-1 - k)
    # Left window for c starts at c-k  -> decr_k[0 : n-2k]
    # Right window for c starts at c   -> incr_k[k : k + (n-2k)]
    num_cands = n - 2 * k
    left_ok = decr_k[0:num_cands]
    right_ok = incr_k[k : k + num_cands]
    mask = left_ok & right_ok

    # Map mask positions (0..num_cands-1) back to r indices c = i + k
    return np.nonzero(mask)[0] + k


# ---------------------------
# Worker: one DV for fixed (et0, coe_op[si])
# ---------------------------
def _dv_worker(et0, coe_op_si, dv_mag, prop_days, max_alt, debug_pc):
    """
    Returns a dict with results for this DV or None if rejected/failed.
    """
    pnt.set_lupnt_epoch(0.0)

    # Time grid: every 10 minutes for prop_days
    tspan = np.linspace(0, prop_days * 24 * 3600, prop_days * 24 * 6 + 1)
    et = et0 - tspan  # dense output times (backward array as in your original)

    # Initial lunar orbit state -> MCI at et0
    rv0_op = pnt.classical_to_cart(coe_op_si, pnt.GM_MOON)
    rv0_mci = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)

    # LOI burn opposite flight direction
    vhat = rv0_mci[3:6] / np.linalg.norm(rv0_mci[3:6])
    dv_vec = -dv_mag * vhat
    rv0_post = rv0_mci - np.hstack((np.zeros(3), dv_vec))

    # Propagate
    xprop_eci, xprop_mci, terminate_fail = propagate_patched_conics(et, rv0_post, debug_pc)

    if terminate_fail:
        return None

    # Find TLI index
    r_earth = np.linalg.norm(xprop_eci[:, 0:3], axis=1)
    idx_close = np.where(r_earth < (pnt.R_EARTH + 500.0e3))[0]

    if len(idx_close) > 0:
        idx_TLI = int(idx_close[0])
    else:
        min_alt = float(np.min(r_earth) - pnt.R_EARTH)
        if min_alt > max_alt:
            return None
        idx_TLI = int(np.argmin(r_earth))

    # TLI info
    et_TLI = float(et[idx_TLI])
    tof = float(et0 - et_TLI)  # [s]

    rv_TLI = xprop_eci[idx_TLI, :]
    r_earth = np.linalg.norm(xprop_eci[:idx_TLI, 0:3], axis=1)
    r_moon = np.linalg.norm(xprop_mci[:idx_TLI, 0:3], axis=1)
    r_TLI = float(np.linalg.norm(rv_TLI[0:3]))
    v_TLI = float(np.linalg.norm(rv_TLI[3:6]))

    # count the number of perigee and perlilune (decrease->increase)
    # identify as perigee if r_earth is decreasing before (5 epochs) and increasing after (5 epochs)
    perigee_idx = k_step_minima_idxs(r_earth, k=5)
    perlilune_idx = k_step_minima_idxs(r_moon, k=5)

    # count as perigee only if altitude < 10000 km
    perigee_idx = perigee_idx[r_earth[perigee_idx] < 1e8]
    perlilune_idx = perlilune_idx[r_moon[perlilune_idx] < 0.5e8]
    n_perigee = len(perigee_idx)
    n_perlilune = len(perlilune_idx)

    # convert to ECEF frame
    rv_TLI_ecef = pnt.convert_frame(et_TLI, rv_TLI, pnt.ECI, pnt.ECEF, rotate_only=True)
    coe = pnt.cart_to_classical(rv_TLI_ecef, pnt.GM_EARTH)

    alt_TLI = (r_TLI - pnt.R_EARTH) / 1e3  # km
    i_TLI = float(coe[2] * 180.0 / np.pi)
    C3_TLI = float(v_TLI**2 - 2.0 * pnt.GM_EARTH / r_TLI)
    C3_TLI = C3_TLI / 1e6  # m^2/s^2 -> km^2/s^2

    # Expand coe_op_si into columns to keep CSV clean
    coe_cols = list(map(float, coe_op_si[:6])) + (
        [float(coe_op_si[6])] if len(coe_op_si) > 6 else []
    )
    row = {
        "et0": float(et0),
        "et_TLI": et_TLI,
        "tof_day": tof / 86400.0,
        "dv_ms": float(dv_mag),
        "C3_mk2s2": C3_TLI,
        "alt_TLI_km": float(alt_TLI),
        "i_TLI_deg": i_TLI,
        "max_r_earth_km": float(np.max(r_earth) / 1000),
        "n_perigee": n_perigee,
        "n_perlilune": n_perlilune,
        "coe_a_km": coe_cols[0] / 1e3,
        "coe_e": coe_cols[1],
        "coe_i_rad": coe_cols[2],
        "coe_Omega_rad": coe_cols[3],
        "coe_omega_rad": coe_cols[4],
        "coe_f_rad": coe_cols[5] if len(coe_cols) > 5 else 0.0,
    }
    return row


# ---------------------------
# Public API: parallel over DVs only; one file per (et0, sat)
# ---------------------------
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import csv
import numpy as np
import sys


def search_low_energy_transfer(
    et0s,
    coes_op,
    dvvec,
    out_dir,
    prop_days=160,
    max_alt=1e5,
    debug=False,
    debug_pc=False,
    n_jobs=None,
):
    """
    Parallelizes only over dvvec for each (et0, sat).
    Saves one CSV file per (et0, sat) into out_dir.
    Prints progress as each DV finishes.
    """
    coes_op = np.array(coes_op, copy=True)
    if coes_op.shape[1] >= 6:
        coes_op[:, 5] = 0.0  # force perilune

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_sat = coes_op.shape[0]
    written = []

    for et0 in et0s:
        if debug:
            print(f"[et0={et0}] processing {n_sat} satellites ...")

        for si in range(n_sat):
            coe_op_si = coes_op[si, :]
            total = len(dvvec)
            done = 0
            rows = []
            fname = out_dir / f"lowE_{et0:.2f}_sat{si:02d}.csv"

            if os.path.exists(fname):
                print(f"  sat {si:02d}: {fname} exists, skipping.")
                continue
            else:
                # --- launch DV workers ---
                with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                    futures = [
                        ex.submit(
                            _dv_worker,  # same worker as in previous example
                            et0,
                            coe_op_si,
                            float(dv_mag),
                            prop_days,
                            max_alt,
                            debug_pc,
                        )
                        for dv_mag in dvvec
                    ]

                    for fut in as_completed(futures):
                        res = fut.result()
                        done += 1
                        # print incremental progress
                        pct = 100.0 * done / total
                        sys.stdout.write(
                            f"\r  sat {si:02d}: {done}/{total} DV cases done ({pct:4.1f}%)"
                        )
                        sys.stdout.flush()
                        if res is not None:
                            rows.append(res)

                sys.stdout.write("\n")  # newline after this sat
                sys.stdout.flush()

                if not rows:
                    # create a empty file to mark no successful transfers
                    with open(fname, "w") as f:
                        rows_keys = [
                            "et0",
                            "et_TLI",
                            "tof_day",
                            "dv_ms",
                            "C3_mk2s2",
                            "alt_TLI_km",
                            "i_TLI_deg",
                            "max_r_earth_km",
                            "n_perigee",
                            "n_perlilune",
                            "coe_a_km",
                            "coe_e",
                            "coe_i_rad",
                            "coe_Omega_rad",
                            "coe_omega_rad",
                            "coe_f_rad",
                        ]
                        w = csv.DictWriter(f, fieldnames=rows_keys)
                        w.writeheader()
                    written.append(str(fname))

                    if debug:
                        print(f"  sat {si:02d}: no successful transfers.")
                    continue
                else:
                    rows.sort(key=lambda r: r["dv_ms"])
                    with open(fname, "w", newline="") as f:
                        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                        w.writeheader()
                        w.writerows(rows)
                    written.append(str(fname))

                if debug:
                    print(f"  sat {si:02d}: wrote {fname} ({len(rows)} rows)")

    if debug:
        print(f"Done. Wrote {len(written)} files to {out_dir}")
    return written


def get_transfer_trajectory(et0, coe_op_si, dv_mag, prop_days=160, debug_pc=False):
    """
    Returns a dict with results for this DV or None if rejected/failed.

    Parameters
    ----------
    et0 : float
        Epoch time of lunar orbit insertion (UTC seconds past J2000).
    coe_op_si : array_like [6] or [7]
        Classical orbital elements of the initial lunar orbit at et0 in Moon-centered
        perifocal frame (a [m], e, i [rad], Omega [rad], omega [rad], f [rad])
    dv_mag : float
        Magnitude of the LOI delta-V in m/s (applied opposite to velocity vector).
    prop_days : int
        Number of days to propagate backward from et0.
    """
    pnt.set_lupnt_epoch(0.0)

    # Time grid: every 10 minutes for prop_days
    tspan = np.linspace(0, prop_days * 24 * 3600, prop_days * 24 * 6 + 1)
    et = et0 - tspan  # dense output times (backward array as in your original)

    # Initial lunar orbit state -> MCI at et0
    rv0_op = pnt.classical_to_cart(coe_op_si, pnt.GM_MOON)
    rv0_mci = pnt.convert_frame(et0, rv0_op, pnt.MOON_OP, pnt.MOON_CI)

    # LOI burn opposite flight direction
    vhat = rv0_mci[3:6] / np.linalg.norm(rv0_mci[3:6])
    dv_vec = -dv_mag * vhat
    rv0_post = rv0_mci - np.hstack((np.zeros(3), dv_vec))

    # Propagate
    xprop_eci, xprop_mci, terminate_fail = propagate_patched_conics(et, rv0_post, debug_pc)
    xprop_serot = eci_to_serot(et, xprop_eci)

    # Find TLI index
    r_earth = np.linalg.norm(xprop_eci[:, 0:3], axis=1)
    r_moon = np.linalg.norm(xprop_mci[:, 0:3], axis=1)
    idx_close = np.where(r_earth < (pnt.R_EARTH + 500.0e3))[0]

    if len(idx_close) > 0:
        idx_TLI = int(idx_close[0])
    else:
        min_alt = float(np.min(r_earth) - pnt.R_EARTH)
        idx_TLI = int(np.argmin(r_earth))

    xprop_eci = xprop_eci[: idx_TLI + 1, :]
    xprop_mci = xprop_mci[: idx_TLI + 1, :]
    xprop_serot = xprop_serot[: idx_TLI + 1, :]

    return xprop_eci, xprop_mci, xprop_serot
