import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Dict, Tuple, Optional, List
import pylupnt as pnt

Vec3 = np.ndarray


# -------------------------------
# Required user-provided adapters
# -------------------------------
# You must provide these three functions (thin wrappers around your tools):
def rtn_to_eci_mat(x_eci: np.ndarray):
    rotmat = np.zeros((3, 3))

    return rotmat


def next_event(t0: float, state0: np.ndarray, event: str) -> Tuple[float, np.ndarray]:
    """Return (t_event, state_event) for event in {"apolune","node"} starting from (t0,state0)."""


def propagate_to(
    t0: float,
    state0: np.ndarray,
    t1: float,
    dv_rtn_at_t0: Optional[np.ndarray] = None,
    dyn: pnt.NBodyDynamics = None,
) -> np.ndarray:
    """Propagate from (t0,state0) to t1. If dv_rtn_at_t0 is not None, apply that RTN Δv at t0."""
    if dv_rtn_at_t0 is not None:
        r2i_mat = rtn_to_eci_mat(state0)
        dv_i = r2i_mat @ dv_rtn_at_t0
        dstate = np.zeros([0, 0, 0, dv_i[0], dv_i[1], dv_i[2]])
        state0 = state0 + dstate
    state_prop = dyn.propagate(state0, t0, t1)

    return state_prop


def osculating_elements(state: np.ndarray) -> Dict[str, float]:
    """Return {"hp": km, "period": s, "i": rad, "RAAN": rad, ...} from the state."""

    coe = pnt.classical_to_cart(state, pnt.GM_MOON)
    a = coe[0]
    e = coe[1]
    inc = coe[2]
    raan = coe[3]
    rp = a * (1 - e)
    hp = rp - pnt.R_MOON
    period = 2 * np.pi * np.sqrt(pnt.GM_MOON / (a**3))

    osc_dict = {
        "hp": hp,
        "period": period,
        "i": inc,
        "raan": raan,
    }

    return osc_dict


# -------------------------------
# Utilities
# -------------------------------
def norm(v: Vec3) -> float:
    return float(np.linalg.norm(v))


def lstsq_weighted(S: np.ndarray, y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Solve min ||W(S dv - y)||_2, with W=diag(w)."""
    W = np.diag(w)
    return np.linalg.lstsq(W @ S, W @ y, rcond=None)[0]


def wrap_pm_pi(angle: float) -> float:
    a = (angle + np.pi) % (2 * np.pi) - np.pi
    return a


@dataclass
class SKTargets:
    # targets
    hp_km: Tuple[float, float]  # (min, max)
    period_s: Optional[float] = None  # desired period, if used
    slot_angle_rad: Optional[float] = None  # per-sat desired mean-longitude offset vs plane ref
    inc_rad: Optional[float] = None
    raan_ref_rad: Optional[float] = None  # use differences vs plane reference
    # deadbands (tolerances)
    db_hp_km: float = 10.0
    db_period_s: float = 20.0
    db_slot_rad: float = np.deg2rad(5.0)
    db_inc_rad: float = np.deg2rad(0.02)
    db_raan_rad: float = np.deg2rad(0.05)


@dataclass
class SKWeights:
    w_hp: float = 5.0
    w_period: float = 3.0
    w_slot: float = 1.0
    w_inc: float = 2.0
    w_raan: float = 2.0


@dataclass
class SKConfig:
    dv_bump_mps: float = 0.01  # finite-diff test burn (m/s)
    min_apolune_interval_days: float = 7.0
    min_node_interval_days: float = 7.0
    combine_inplane: bool = True  # do one combined in-plane burn at apolune
    combine_cross: bool = True  # one normal burn near node
    max_iter_per_event: int = 1  # usually 1 is enough for small errors
    rng: np.random.Generator = field(default_factory=lambda: np.random.default_rng(1))


@dataclass
class Burn:
    t: float
    event: str  # "apolune" or "node"
    dv_rtn_mps: Vec3  # RTN components [dr, dt, dn]
    dv_mag_mps: float


@dataclass
class SKLog:
    burns: List[Burn] = field(default_factory=list)
    dv_total_mps: float = 0.0


# -------------------------------
# Sensitivity via finite-difference
# -------------------------------
def build_sensitivity_inplane(
    t0: float,
    state0: np.ndarray,
    next_event: Callable,
    propagate_to: Callable,
    osculating_elements: Callable,
    targets: SKTargets,
    use_slot: bool,
    dv_bump: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build S (rows = [hp, period, slot?], cols = [r,t,n]) at APO,
    and the current error vector y.
    """
    # Go to apolune without burn
    t_ap, st_ap = next_event(t0, state0, "apolune")
    el = osculating_elements(st_ap)

    rows = []
    y = []

    # hp error: clamp inside band -> aim for center
    hp_c = 0.5 * (targets.hp_km[0] + targets.hp_km[1])
    rows.append(lambda e: e["hp"])
    y.append(hp_c - el["hp"])

    # period error
    if targets.period_s is not None:
        rows.append(lambda e: e["period"])
        y.append(targets.period_s - el["period"])

    # slot (mean-longitude offset)
    if use_slot and (targets.slot_angle_rad is not None) and ("mean_longitude" in el):
        # expect user to supply plane reference separately; here we target the provided slot directly
        rows.append(lambda e: e["mean_longitude"])
        # y will be filled outside with provided slot target vs current
        # (we leave as current since absolute targeting requires reference; user can inject externally)
        y.append(0.0)

    m = len(rows)
    S = np.zeros((m, 3))

    # Bump along r,t,n at apolune
    for k, axis in enumerate(np.eye(3)):
        dv = dv_bump * axis
        st_plus = propagate_to(t_ap, st_ap, t_ap, dv_rtn_at_t0=dv)  # instantaneous at apolune
        # propagate to next apolune to capture effect on mean quantities
        t_ap2, st_ap2 = next_event(t_ap, st_plus, "apolune")
        elp = osculating_elements(st_ap2)
        # likewise nominal no-burn to same endpoint
        t_ap2n, st_ap2n = next_event(t_ap, st_ap, "apolune")
        eln = osculating_elements(st_ap2n)

        for i, f in enumerate(rows):
            S[i, k] = (f(elp) - f(eln)) / dv_bump

    return S, np.array(y)


def build_sensitivity_cross(
    t0: float,
    state0: np.ndarray,
    next_event: Callable,
    propagate_to: Callable,
    osculating_elements: Callable,
    targets: SKTargets,
    dv_bump: float,
    dynamics: pnt.NBodyDynamics,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build S (rows = [inc, raan], cols = [r,t,n]) near NODE (ascending).
    """
    t_nd, st_nd = next_event(t0, state0, "node")
    el = osculating_elements(st_nd)

    rows = [lambda e: e["i"], lambda e: e["RAAN"]]
    y = []
    if targets.inc_rad is not None:
        y.append(targets.inc_rad - el["i"])
    else:
        y.append(0.0)

    if targets.raan_ref_rad is not None:
        dRAAN = wrap_pm_pi(targets.raan_ref_rad - el["RAAN"])
        y.append(dRAAN)
    else:
        y.append(0.0)

    m = len(rows)
    S = np.zeros((m, 3))

    for k, axis in enumerate(np.eye(3)):
        dv = dv_bump * axis
        st_plus = propagate_to(t_nd, st_nd, t_nd, dv_rtn_at_t0=dv)
        # propagate to next node to evaluate mean plane change
        t_nd2, st_nd2 = next_event(t_nd, st_plus, "node")
        elp = osculating_elements(st_nd2)

        t_nd2n, st_nd2n = next_event(t_nd, st_nd, "node")
        eln = osculating_elements(st_nd2n)

        for i, f in enumerate(rows):
            d = f(elp) - f(eln)
            if i == 1:  # RAAN wrap
                d = wrap_pm_pi(d)
            S[i, k] = d / dv_bump

    return S, np.array(y)


# -------------------------------
# Controller
# -------------------------------
class StationKeeper:
    def __init__(
        self, targets: SKTargets, weights: SKWeights, cfg: SKConfig, dynamics: pnt.NBodyDynamics
    ):
        self.targets = targets
        self.weights = weights
        self.cfg = cfg
        self.dynamics = dynamics

    def plan_apolune_burn(
        self, t0, state0, next_event, propagate_to, osculating_elements, slot_error_rad: float = 0.0
    ) -> Optional[Burn]:
        """Compute combined in-plane Δv at apolune to fix hp/period/(slot)."""
        S, y = build_sensitivity_inplane(
            t0,
            state0,
            next_event,
            propagate_to,
            osculating_elements,
            self.targets,
            use_slot=(abs(slot_error_rad) > 0),
            dv_bump=self.cfg.dv_bump_mps,
        )

        # fill slot target if requested (append weight if present)
        w = []
        # hp
        w.append(self.weights.w_hp)
        # period
        if self.targets.period_s is not None:
            w.append(self.weights.w_period)
        # slot
        if abs(slot_error_rad) > 0:
            # Overwrite last y element with desired slot correction
            y[-1] = slot_error_rad
            w.append(self.weights.w_slot)

        dv_rtn = lstsq_weighted(S, y, np.array(w))
        if np.allclose(dv_rtn, 0.0, atol=1e-6):
            return None

        t_ap, _ = next_event(t0, state0, "apolune")
        return Burn(t=t_ap, event="apolune", dv_rtn_mps=dv_rtn, dv_mag_mps=norm(dv_rtn))

    def plan_node_burn(
        self, t0, state0, next_event, propagate_to, osculating_elements
    ) -> Optional[Burn]:
        """Compute cross-track Δv at node to fix i/RAAN (mostly normal axis)."""
        S, y = build_sensitivity_cross(
            t0,
            state0,
            next_event,
            propagate_to,
            osculating_elements,
            self.targets,
            dv_bump=self.cfg.dv_bump_mps,
        )
        # weights: [i, RAAN]
        w = np.array([self.weights.w_inc, self.weights.w_raan])
        dv_rtn = lstsq_weighted(S, y, w)
        if np.allclose(dv_rtn, 0.0, atol=1e-6):
            return None
        t_nd, _ = next_event(t0, state0, "node")
        return Burn(t=t_nd, event="node", dv_rtn_mps=dv_rtn, dv_mag_mps=norm(dv_rtn))

    def run(
        self,
        t_start: float,
        t_end: float,
        state0: np.ndarray,
        next_event: Callable,
        propagate_to: Callable,
        osculating_elements: Callable,
        plane_ref: Optional[Dict[str, float]] = None,  # e.g., {"mean_longitude": ... , "RAAN": ...}
        last_apolune_burn_t: float = -np.inf,
        last_node_burn_t: float = -np.inf,
    ) -> Tuple[np.ndarray, SKLog]:
        """
        Event-driven loop: checks weekly for node trim; each apolune for in-plane trim.
        Returns final state and SK log.
        """
        log = SKLog()
        t = t_start
        state = state0.copy()

        while t < t_end:
            # 1) Apolune opportunity
            t_ap, st_ap = next_event(t, state, "apolune")
            if t_ap > t_end:  # nothing else to do
                break

            el_ap = osculating_elements(st_ap)
            # Decide if in-plane correction is needed
            needs_hp = not (self.targets.hp_km[0] <= el_ap["hp"] <= self.targets.hp_km[1])
            needs_T = (self.targets.period_s is not None) and (
                abs(el_ap["period"] - self.targets.period_s) > self.targets.db_period_s
            )
            needs_slot = False
            slot_err = 0.0
            if (
                (self.targets.slot_angle_rad is not None)
                and ("mean_longitude" in el_ap)
                and (plane_ref is not None)
            ):
                slot_err = wrap_pm_pi(
                    (plane_ref["mean_longitude"] + self.targets.slot_angle_rad)
                    - el_ap["mean_longitude"]
                )
                needs_slot = abs(slot_err) > self.targets.db_slot_rad

            do_inplane = (needs_hp or needs_T or needs_slot) and (
                t_ap - last_apolune_burn_t >= self.cfg.min_apolune_interval_days * 86400.0
            )

            if do_inplane:
                burn = self.plan_apolune_burn(
                    t,
                    state,
                    next_event,
                    propagate_to,
                    osculating_elements,
                    slot_error_rad=slot_err if needs_slot else 0.0,
                )
                if burn and burn.dv_mag_mps > 1e-5:
                    # execute burn at t_ap
                    state = propagate_to(t, state, t_ap)  # coast to apolune
                    state = propagate_to(t_ap, state, t_ap, dv_rtn_at_t0=burn.dv_rtn_mps)
                    last_apolune_burn_t = t_ap
                    log.burns.append(burn)
                    log.dv_total_mps += burn.dv_mag_mps
                    t = t_ap
                else:
                    # no burn; just move on to t_ap
                    state = propagate_to(t, state, t_ap)
                    t = t_ap
            else:
                # Just move to apolune time
                state = propagate_to(t, state, t_ap)
                t = t_ap

            # 2) Weekly node opportunity (check from current time)
            t_nd, st_nd = next_event(t, state, "node")
            if t_nd <= t_end and (
                t_nd - last_node_burn_t >= self.cfg.min_node_interval_days * 86400.0
            ):
                el_nd = osculating_elements(st_nd)
                needs_i = (self.targets.inc_rad is not None) and (
                    abs(el_nd["i"] - self.targets.inc_rad) > self.targets.db_inc_rad
                )
                needs_O = (self.targets.raan_ref_rad is not None) and (
                    abs(wrap_pm_pi(el_nd["RAAN"] - self.targets.raan_ref_rad))
                    > self.targets.db_raan_rad
                )

                if needs_i or needs_O:
                    burn = self.plan_node_burn(
                        t, state, next_event, propagate_to, osculating_elements
                    )
                    if burn and burn.dv_mag_mps > 1e-5:
                        state = propagate_to(t, state, t_nd)  # coast
                        state = propagate_to(t_nd, state, t_nd, dv_rtn_at_t0=burn.dv_rtn_mps)
                        last_node_burn_t = t_nd
                        log.burns.append(burn)
                        log.dv_total_mps += burn.dv_mag_mps
                        t = t_nd
                        continue  # next loop from node time

            # If we didn’t burn at node (or node is far), advance a fraction of an orbit to avoid infinite loop
            # Here we step to just after apolune to look for the next apolune
            t = t + 0.1 * el_ap["period"]

        # Final coast to t_end if needed
        if t < t_end:
            state = propagate_to(t, state, t_end)

        return state, log
