# --- blt_surrogate.py ---
import numpy as np
import spiceypy as sp
from datetime import datetime, timezone


def sun_phase_sem(r_earth_to_moon_km, r_earth_to_sun_km):
    """
    Sun–Earth–Moon phase angle φ in [0, 2π), measured between Earth→Moon and Earth→Sun.
    """
    uM = r_earth_to_moon_km / (np.linalg.norm(r_earth_to_moon_km) + 1e-16)
    uS = r_earth_to_sun_km / (np.linalg.norm(r_earth_to_sun_km) + 1e-16)
    c = np.clip(np.dot(uM, uS), -1.0, 1.0)
    s = np.linalg.norm(np.cross(uM, uS))
    phi = float(np.arctan2(s, c)) % (2 * np.pi)
    return phi


def dv_loi_blt_surrogate(
    dv_loi_direct_kms, tof_days, sun_phase_rad, savings_table=None, dv_corr_ms=50.0
):
    """
    Replace direct LOI ΔV with a BLT-style surrogate:
      dv_LOI_BLT = max(0, dv_LOI_direct - savings(TOF, sun_phase)) + dv_corr
    Defaults give ~0.2–0.6 km/s savings growing with TOF (60–120 d), plus small corrections.
    """
    if savings_table is None:
        # heuristic: more TOF → more savings (saturates ~0.6 km/s), ±0.05 with phase
        base = 0.20
        extra = 0.40 * np.clip((tof_days - 60.0) / 60.0, 0.0, 1.0)  # 60→120 d
        phase_mod = 0.05 * np.cos(sun_phase_rad)
        savings = np.clip(base + extra + phase_mod, 0.15, 0.60)
    else:
        savings = float(savings_table(tof_days, sun_phase_rad))

    dv_corr = dv_corr_ms / 1000.0
    return max(0.0, dv_loi_direct_kms - savings) + dv_corr


# --- spice_adapter.py (extend your existing file) ---


def load_spice_kernels(kernel_paths):
    for k in kernel_paths:
        sp.furnsh(k)


def _dt_to_et(t_dt: datetime) -> float:
    if t_dt.tzinfo is None:
        t_dt = t_dt.replace(tzinfo=timezone.utc)
    return sp.utc2et(t_dt.strftime("%Y-%m-%dT%H:%M:%S"))


def moon_state_spice(t_dt: datetime):
    et = _dt_to_et(t_dt)
    state, _ = sp.spkezr("MOON", et, "J2000", "NONE", "EARTH")
    return np.asarray(state[:3]), np.asarray(state[3:])


def sun_state_spice(t_dt: datetime):
    et = _dt_to_et(t_dt)
    state, _ = sp.spkezr("SUN", et, "J2000", "NONE", "EARTH")
    return np.asarray(state[:3]), np.asarray(state[3:])


# --- dv_walker_blt.py ---
import numpy as np
from datetime import timedelta
from typing import Callable, Tuple, Optional, List, Dict

from blt_surrogate import dv_loi_blt_surrogate, sun_phase_sem
from dv_walker_elfo import (
    ELFO,
    dv_earth_to_elfo_elts,
    ephem_knobs_from_moon_state,
    R_E,
    R_M,
    A_MOON,
    MU_E,
)


def dv_for_walker_constellation_blt(
    h_leo_km: float,
    elfo: ELFO,
    t0_utc,
    tof_days: float,
    targets_OM_list: List[Tuple[float, float, int, int]],
    moon_state_func: Callable = None,  # e.g., spice_adapter.moon_state_spice
    sun_state_func: Callable = None,  # e.g., spice_adapter.sun_state_spice
    use_free_phasing: bool = True,
    N_revs_phase: int = 1,
    delta_i_E_rad: float = 0.0,
    delta_argp_rad: float = 0.0,
    use_blt: bool = True,
    blt_savings_table: Optional[Callable] = None,
    blt_corr_ms: float = 50.0,
) -> Tuple[float, List[Dict]]:
    """
    Same as your previous batch driver, but with an optional BLT/WSB LOI surrogate.
    """
    assert (
        moon_state_func is not None and sun_state_func is not None
    ), "Pass Moon and Sun ephemeris functions."

    t_arrival = t0_utc + timedelta(days=float(tof_days))

    # Ephemerides at arrival
    rM_eci, vM_eci = moon_state_func(t_arrival)
    rS_eci, vS_eci = sun_state_func(t_arrival)

    # Arrival geometry knobs (same as before)
    r_mag = float(np.linalg.norm(rM_eci))
    a_T = 0.5 * ((R_E + h_leo_km) + A_MOON)
    vT_at_rM = float(np.sqrt(MU_E * (2.0 / r_mag - 1.0 / a_T)))
    cos_phase, i_arrival, _vinf, raan_arrival = ephem_knobs_from_moon_state(
        rM_eci, vM_eci, vT_at_rM
    )

    # Sun–Earth–Moon phase (for BLT savings surrogate)
    phi_sun = sun_phase_sem(rM_eci, rS_eci)

    # ELFO mean motion for phase-time mapping
    rpM, raM, aL, eL = elfo.rp_ra()
    n_elfo = np.sqrt(4902.800066 / aL**3)

    dv_max = 0.0
    details = []

    for raan_tgt, M_tgt, p_idx, s_idx in targets_OM_list:
        # ΔΩ at apolune
        dO = (raan_tgt - raan_arrival + np.pi) % (2 * np.pi) - np.pi

        # Core (direct) DV and parts, including dv_LOI_direct
        dv_direct, parts = dv_earth_to_elfo_elts(
            h_leo_km=h_leo_km,
            elfo=elfo,
            cos_phase=cos_phase,
            i_arrival_rad=i_arrival,
            delta_raan_rad=dO,
            delta_argp_rad=delta_argp_rad,
            delta_i_E_rad=delta_i_E_rad,
        )
        dv_LOI_direct = parts["dv_LOI"]

        # BLT swap for LOI
        if use_blt:
            dv_LOI_blt = dv_loi_blt_surrogate(
                dv_loi_direct_kms=dv_LOI_direct,
                tof_days=float(tof_days),
                sun_phase_rad=phi_sun,
                savings_table=blt_savings_table,
                dv_corr_ms=blt_corr_ms,
            )
            dv_core = dv_direct - dv_LOI_direct + dv_LOI_blt
        else:
            dv_core = dv_direct
            dv_LOI_blt = dv_LOI_direct

        # Phase handling (same as before)
        if use_free_phasing:
            delta_t_sec = M_tgt / n_elfo
            dv_phase = 0.0
            onstation_epoch = t_arrival + timedelta(seconds=float(delta_t_sec))
        else:
            # reuse your two-impulse apolune phasing function if you want fixed-epoch
            from dv_walker_elfo import dv_phase_two_impulses_apolune

            dv_phase, _ = dv_phase_two_impulses_apolune(
                aL, rpM, raM, delta_M_rad=M_tgt, N_revs=N_revs_phase
            )
            onstation_epoch = t_arrival

        dv_total = dv_core + dv_phase
        dv_max = max(dv_max, dv_total)

        details.append(
            dict(
                plane=p_idx,
                sat=s_idx,
                raan_target=raan_tgt,
                M_target=M_tgt,
                dRAAN=dO,
                dv_total=dv_total,
                dv_core=dv_core,
                dv_phase=dv_phase,
                dv_LOI_direct=dv_LOI_direct,
                dv_LOI_blt=dv_LOI_blt,
                sun_phase_rad=phi_sun,
                arrival_epoch=onstation_epoch,
                **parts
            )
        )

    return dv_max, details
