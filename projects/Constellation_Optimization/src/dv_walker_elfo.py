import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple, Optional, List, Dict
from datetime import datetime, timedelta
import pylupnt as pnt

# ---------- constants ----------
MU_E = pnt.GM_EARTH  # km^3/s^2
MU_M = pnt.GM_MOON
R_E = pnt.R_EARTH  # km
R_M = pnt.R_MOON  # km
A_MOON = 384400.0  # km
N_MOON = np.sqrt(MU_E / A_MOON**3)
V_MOON = np.sqrt(MU_E / A_MOON)
J2000_Z = np.array([0.0, 0.0, 1.0])


# ---------- utils ----------
def _norm(v):
    return float(np.linalg.norm(v))


def _uhat(v):
    return v / (np.linalg.norm(v) + 1e-16)


def _prograde_tangent(r_hat, plane_normal_hat, v_ref):
    t_hat = _uhat(np.cross(plane_normal_hat, r_hat))
    return t_hat if np.dot(t_hat, v_ref) >= 0 else -t_hat


@dataclass
class ELFO:
    a_km: Optional[float] = None
    e: Optional[float] = None
    rp_km: Optional[float] = None
    ra_km: Optional[float] = None
    i_rad: float = 0.0

    def rp_ra(self):
        if self.a_km is not None and self.e is not None:
            rp = self.a_km * (1 - self.e)
            ra = self.a_km * (1 + self.e)
            return rp, ra, self.a_km, self.e
        if self.rp_km is not None and self.ra_km is not None:
            a = 0.5 * (self.rp_km + self.ra_km)
            e = (self.ra_km - self.rp_km) / (self.ra_km + self.rp_km + 1e-30)
            return self.rp_km, self.ra_km, a, e
        raise ValueError("Provide (a,e) or (rp,ra).")


# very light circular Moon ephemeris (swap for SPICE in real use)
def circular_moon_state(t: datetime, t0=datetime(2000, 1, 1, 12, 0, 0)):
    dt = (t - t0).total_seconds()
    th = N_MOON * dt
    c, s = np.cos(th), np.sin(th)
    r = np.array([A_MOON * c, A_MOON * s, 0.0])
    v = np.array([-A_MOON * N_MOON * s, A_MOON * N_MOON * c, 0.0])
    return r, v


def ephem_knobs_from_moon_state(rM, vM, vT_mag_at_rM, zL_hat=J2000_Z):
    r_hat = _uhat(rM)
    plane_n = _uhat(np.cross(rM, vM))
    t_hat_T = _prograde_tangent(r_hat, plane_n, vM)
    vT_vec = vT_mag_at_rM * t_hat_T
    vinf_vec = vT_vec - vM
    cos_phase = float(np.clip(np.dot(_uhat(vT_vec), _uhat(vM)), -1, 1))
    h_inf = np.cross(rM, vinf_vec)
    i_arrival = float(np.arccos(np.clip(np.dot(_uhat(h_inf), _uhat(zL_hat)), -1, 1)))
    # RAAN of arrival plane wrt J2000 equator (for ΔΩ bookkeeping)
    n_vec = _uhat(h_inf)
    node_line = np.cross(J2000_Z, n_vec)
    if _norm(node_line) < 1e-12:
        raan_arrival = 0.0
    else:
        node_hat = _uhat(node_line)
        raan_arrival = float(np.arctan2(node_hat[1], node_hat[0])) % (2 * np.pi)
    return cos_phase, i_arrival, vinf_vec, raan_arrival


# ---------- single-sat DV core ----------
def dv_earth_to_elfo_elts(
    h_leo_km: float,
    elfo: ELFO,
    cos_phase: float,
    i_arrival_rad: float,
    delta_raan_rad: float = 0.0,
    delta_argp_rad: float = 0.0,
    delta_i_E_rad: float = 0.0,
) -> Tuple[float, Dict[str, float]]:
    rpM, raM, aL, eL = elfo.rp_ra()
    r_p = R_E + h_leo_km
    v_LEO = np.sqrt(MU_E / r_p)
    a_T = 0.5 * (r_p + A_MOON)
    v_pT = np.sqrt(MU_E * (2 / r_p - 1 / a_T))
    dv_TLI = v_pT - v_LEO

    v_T_at_M = np.sqrt(MU_E * (2 / A_MOON - 1 / a_T))
    v_inf_sq = v_T_at_M**2 + V_MOON**2 - 2 * v_T_at_M * V_MOON * cos_phase
    v_inf = float(np.sqrt(max(v_inf_sq, 0.0)))

    v_hyp_p = np.sqrt(v_inf**2 + 2 * MU_M / rpM)
    v_L_p = np.sqrt(MU_M * (2 / rpM - 1 / aL))
    dv_LOI = abs(v_hyp_p - v_L_p)

    v_L_a = np.sqrt(MU_M * (2 / raM - 1 / aL))
    delta_i_M = abs(elfo.i_rad - i_arrival_rad)
    theta_plane = np.hypot(delta_i_M, delta_raan_rad)
    dv_plane_M = 2 * v_L_a * np.sin(0.5 * theta_plane)
    dv_argp_M = 2 * v_L_a * np.sin(0.5 * abs(delta_argp_rad))
    dv_plane_E = 2 * v_LEO * np.sin(0.5 * delta_i_E_rad)

    dv_total = dv_TLI + dv_LOI + dv_plane_E + dv_plane_M + dv_argp_M
    parts = dict(
        dv_TLI=dv_TLI,
        dv_LOI=dv_LOI,
        dv_plane_E=dv_plane_E,
        dv_plane_M=dv_plane_M,
        dv_argp_M=dv_argp_M,
        v_inf=v_inf,
        rpM=rpM,
        raM=raM,
        aL=aL,
        eL=eL,
        v_L_a=v_L_a,
    )
    return dv_total, parts


# ---------- apolune 2-impulse phasing to realize ΔM at fixed epoch ----------
def dv_phase_two_impulses_apolune(aL, rpM, raM, delta_M_rad, N_revs=1):
    mu = MU_M
    n_nom = np.sqrt(mu / aL**3)
    T_nom = 2 * np.pi / n_nom
    n_req = n_nom + delta_M_rad / (N_revs * T_nom)
    a_req = (mu / n_req**2) ** (1 / 3)
    r_a = raM
    v_nom_a = np.sqrt(mu * (2 / r_a - 1 / aL))
    v_req_a = np.sqrt(mu * (2 / r_a - 1 / a_req))
    dv = abs(v_req_a - v_nom_a)
    return 2 * dv, N_revs * (2 * np.pi / n_req)  # (dv_total, phasing time)


# ---------- Walker targets (Ω, M) ----------
def walker_delta_targets(P, S, F, raan0_rad=0.0, M0_rad=0.0):
    T = P * S
    out = []
    for p in range(P):
        raan = (raan0_rad + 2 * np.pi * p / P) % (2 * np.pi)
        for s in range(S):
            M = (M0_rad + 2 * np.pi * (s / S + (F * p) / T)) % (2 * np.pi)
            out.append((raan, M, p, s))
    return out


# ---------- batch driver ----------
def dv_for_walker_constellation(
    h_leo_km: float,
    elfo: ELFO,
    t0_utc: datetime,
    tof_days: float,
    targets_OM_list: List[Tuple[float, float, int, int]],
    moon_state_func: Callable[[datetime], Tuple[np.ndarray, np.ndarray]] = circular_moon_state,
    lunar_pole_hat_eci: np.ndarray = J2000_Z,
    use_free_phasing: bool = True,
    fixed_epoch: Optional[datetime] = None,
    N_revs_phase: int = 1,
    delta_i_E_rad: float = 0.0,
    delta_argp_rad: float = 0.0,
) -> Tuple[float, List[Dict]]:
    # common arrival epoch & Moon state
    t_arrival = t0_utc + timedelta(days=float(tof_days))
    rM, vM = moon_state_func(t_arrival)
    r_mag = _norm(rM)
    a_T = 0.5 * (R_E + h_leo_km + A_MOON)
    vT_at_rM = np.sqrt(MU_E * (2 / r_mag - 1 / a_T))
    cos_phase, i_arrival, _vinf, raan_arrival = ephem_knobs_from_moon_state(
        rM, vM, vT_at_rM, zL_hat=lunar_pole_hat_eci
    )

    rpM, raM, aL, eL = elfo.rp_ra()
    n_elfo = np.sqrt(MU_M / aL**3)
    T_elfo = 2 * np.pi / n_elfo

    dv_max = 0.0
    details = []
    for raan_tgt, M_tgt, p_idx, s_idx in targets_OM_list:
        # ΔΩ (wrap to [-π,π])
        dO = (raan_tgt - raan_arrival + np.pi) % (2 * np.pi) - np.pi

        # core DV (Earth->Moon, RAAN at apolune)
        dv_core, parts = dv_earth_to_elfo_elts(
            h_leo_km,
            elfo,
            cos_phase,
            i_arrival,
            delta_raan_rad=dO,
            delta_argp_rad=delta_argp_rad,
            delta_i_E_rad=delta_i_E_rad,
        )

        # phase handling
        if use_free_phasing:
            # realize M via LOI time offset (≈0 DV); report Δt
            delta_t_sec = M_tgt / n_elfo
            dv_phase = 0.0
            onstation_epoch = t_arrival + timedelta(seconds=delta_t_sec)
        else:
            # fixed epoch (everyone same epoch): two-impulse apolune phasing
            dv_phase, _Tphase = dv_phase_two_impulses_apolune(
                aL, rpM, raM, delta_M_rad=M_tgt, N_revs=N_revs_phase
            )
            delta_t_sec = 0.0
            onstation_epoch = fixed_epoch if fixed_epoch else t_arrival

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
                onstation_epoch=onstation_epoch,
                n_elfo=n_elfo,
                T_elfo=T_elfo,
                **parts
            )
        )

    return dv_max, details
