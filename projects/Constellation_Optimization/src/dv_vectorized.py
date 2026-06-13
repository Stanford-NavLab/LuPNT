# ---- dv_vectorized.py ----
import numpy as np
from datetime import datetime, timedelta
from typing import List, Tuple, Optional, Dict, Callable

# import your existing pieces:
from src.dv_walker_elfo import (
    ELFO,
    dv_earth_to_elfo_elts,
    ephem_knobs_from_moon_state,
    dv_phase_two_impulses_apolune,
    walker_delta_targets,
    R_E,
    R_M,
    A_MOON,
    MU_E,
)

# Light re-declarations for standalone snippet:
MU_E = 398600.4418
A_MOON = 384400.0


def _norm(v):
    return float(np.linalg.norm(v))


def _arrival_knobs_for_epoch(
    h_leo_km: float,
    arrival_dt: datetime,
    moon_state_func: Callable[[datetime], Tuple[np.ndarray, np.ndarray]],
):
    rM, vM = moon_state_func(arrival_dt)
    r_mag = _norm(rM)
    a_T = 0.5 * ((6378.1363 + h_leo_km) + A_MOON)
    vT_at_rM = np.sqrt(MU_E * (2.0 / r_mag - 1.0 / a_T))
    cos_phase, i_arrival, _vinf, raan_arrival = ephem_knobs_from_moon_state(rM, vM, vT_at_rM)
    return cos_phase, i_arrival, raan_arrival


def _dv_for_targets_once(
    h_leo_km: float,
    elfo,
    cos_phase: float,
    i_arrival_rad: float,
    raan_arrival: float,
    targets_OM_list: List[Tuple[float, float, int, int]],
    delta_i_E_rad: float,
    delta_argp_rad: float,
):
    # compute core DV for each sat given (cos_phase, i_arrival, ΔΩ)
    dv_core = np.empty(len(targets_OM_list))
    dRAAN = np.empty(len(targets_OM_list))
    parts_list: List[Dict] = []
    for k, (raan_tgt, _M_tgt, _p, _s) in enumerate(targets_OM_list):
        dO = (raan_tgt - raan_arrival + np.pi) % (2 * np.pi) - np.pi
        dRAAN[k] = dO
        dv, parts = dv_earth_to_elfo_elts(
            h_leo_km,
            elfo,
            cos_phase=cos_phase,
            i_arrival_rad=i_arrival_rad,
            delta_raan_rad=dO,
            delta_argp_rad=delta_argp_rad,
            delta_i_E_rad=delta_i_E_rad,
        )
        dv_core[k] = dv
        parts_list.append(parts)
    return dv_core, dRAAN, parts_list


def dv_for_walker_constellation_vec_tof(
    h_leo_km: float,
    elfo,
    t0_utc: datetime,
    tof_days_array: np.ndarray,  # shape (N_toF,)
    targets_OM_list: List[Tuple[float, float, int, int]],
    moon_state_func: Callable[[datetime], Tuple[np.ndarray, np.ndarray]],
    use_free_phasing: bool = True,
    N_revs_phase: int = 1,
    delta_i_E_rad: float = 0.0,
    delta_argp_rad: float = 0.0,
) -> Dict[str, np.ndarray]:
    """
    Vectorized over TOF. Returns dict with:
      'dv_peak' shape (N_toF,), 'dv_core' shape (N_toF, N_sats), 'dRAAN' shape (N_toF, N_sats),
      'arrival_epoch' shape (N_toF,), 'dv_phase' shape (N_toF, N_sats), 'onstation_epoch' shape (N_toF, N_sats)
    """
    N_toF = len(tof_days_array)
    N_sats = len(targets_OM_list)

    dv_core_mat = np.zeros((N_toF, N_sats))
    dRAAN_mat = np.zeros((N_toF, N_sats))
    dv_phase_mat = np.zeros((N_toF, N_sats))
    arrival_epochs = np.empty(N_toF, dtype=object)
    onstation_epochs = np.empty((N_toF, N_sats), dtype=object)

    # mean motion (for phase handling/reporting)
    rpM, raM, aL, eL = elfo.rp_ra()
    n_elfo = np.sqrt(4902.800066 / aL**3)
    T_elfo = 2 * np.pi / n_elfo

    for i, tof in enumerate(tof_days_array):
        arr = t0_utc + timedelta(days=float(tof))
        arrival_epochs[i] = arr

        cos_phase, i_arrival, raan_arrival = _arrival_knobs_for_epoch(
            h_leo_km, arr, moon_state_func
        )

        dv_core, dRAAN, _parts_list = _dv_for_targets_once(
            h_leo_km,
            elfo,
            cos_phase,
            i_arrival,
            raan_arrival,
            targets_OM_list,
            delta_i_E_rad,
            delta_argp_rad,
        )
        dv_core_mat[i, :] = dv_core
        dRAAN_mat[i, :] = dRAAN

        # phase handling per sat
        for k, (_raan_tgt, M_tgt, _p, _s) in enumerate(targets_OM_list):
            if use_free_phasing:
                dt = M_tgt / n_elfo
                dvp = 0.0
                onstation_epochs[i, k] = arr + timedelta(seconds=float(dt))
            else:
                dvp, _ = dv_phase_two_impulses_apolune(
                    aL, rpM, raM, delta_M_rad=M_tgt, N_revs=N_revs_phase
                )
                onstation_epochs[i, k] = arr
            dv_phase_mat[i, k] = dvp

    dv_total_mat = dv_core_mat + dv_phase_mat
    dv_peak = dv_total_mat.max(axis=1)
    return dict(
        dv_peak=dv_peak,
        dv_total=dv_total_mat,
        dv_core=dv_core_mat,
        dv_phase=dv_phase_mat,
        dRAAN=dRAAN_mat,
        arrival_epoch=arrival_epochs,
        onstation_epoch=onstation_epochs,
        T_elfo=T_elfo,
        n_elfo=n_elfo,
    )


def dv_for_walker_constellation_vec_epoch(
    h_leo_km: float,
    elfo,
    t0_list: List[datetime],
    tof_days: float,
    targets_OM_list: List[Tuple[float, float, int, int]],
    moon_state_func: Callable[[datetime], Tuple[np.ndarray, np.ndarray]],
    **kwargs
):
    """
    Vectorized over launch epochs (single TOF).
    Returns the same dictionary shapes but first axis is len(t0_list).
    """
    tof_arr = np.full(len(t0_list), float(tof_days))
    # We evaluate each launch epoch separately (Moon state differs), but use the same vector engine
    outs = []
    for t0 in t0_list:
        outs.append(
            dv_for_walker_constellation_vec_tof(
                h_leo_km,
                elfo,
                t0,
                np.array([tof_arr[0]]),
                targets_OM_list,
                moon_state_func,
                **kwargs
            )
        )
    # Stack results along axis 0
    dv_peak = np.array([o["dv_peak"][0] for o in outs])
    dv_total = np.stack([o["dv_total"][0] for o in outs], axis=0)
    dv_core = np.stack([o["dv_core"][0] for o in outs], axis=0)
    dv_phase = np.stack([o["dv_phase"][0] for o in outs], axis=0)
    dRAAN = np.stack([o["dRAAN"][0] for o in outs], axis=0)
    arrival_epoch = np.array([o["arrival_epoch"][0] for o in outs], dtype=object)
    onstation_epoch = np.stack([o["onstation_epoch"][0] for o in outs], axis=0)
    return dict(
        dv_peak=dv_peak,
        dv_total=dv_total,
        dv_core=dv_core,
        dv_phase=dv_phase,
        dRAAN=dRAAN,
        arrival_epoch=arrival_epoch,
        onstation_epoch=onstation_epoch,
        T_elfo=outs[0]["T_elfo"],
        n_elfo=outs[0]["n_elfo"],
    )
