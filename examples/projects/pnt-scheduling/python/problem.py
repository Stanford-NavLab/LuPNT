import numpy as np
import logging
from dataclasses import dataclass, field
from copy import deepcopy
from typing import ClassVar
from dataclasses import dataclass

# *************************************************************************************************
# Data types
# *************************************************************************************************

ABS_TOL = 1e-6


def reset_id_counters():
    Request._next_id = 0
    ServiceWindow._next_id = 0


class UsrId(int):
    pass


class SatId(int):
    pass


class ReqId(int):
    pass


class WinId(int):
    pass


class Data(float):
    pass


class Energy(float):
    pass


class Time(float):
    pass


@dataclass(frozen=True, repr=True)
class Request:
    usr_id: UsrId  # User id
    ts: float  # Start time
    te: float  # End time
    dur: float  # Duration
    ta: float = 0  # Arrival time
    id: int = field(default_factory=lambda: Request._next_id)  # Request id

    _next_id: ClassVar[int] = 0

    def __post_init__(self):
        Request._next_id += 1


@dataclass(frozen=True, repr=True)
class ServiceWindow:
    sat_id: SatId  # Satellite id
    usr_id: int  # User id
    ts: float  # Start time
    te: float  # End time
    id: WinId = field(default_factory=lambda: ServiceWindow._next_id)  # Window id

    _next_id: ClassVar[WinId] = 0

    def __post_init__(self):
        ServiceWindow._next_id += 1


@dataclass(frozen=True, repr=False)
class State:
    t: list[float]  # Current time
    req: list[Request]  # Last request
    req_times: dict[ReqId, float]  # Request time
    D: list[float]  # Data onboard
    E: list[float]  # Energy onboard

    def __repr__(self):
        s = "State("
        s += f"t={np.round(self.t, 2)}, "
        s += f"D={np.round(self.D, 2)}, "
        s += f"E={np.round(self.E, 2)}, "
        s += f"req={list(req.id if req else None for req in self.req)}, "
        s += f"req_times={list(self.req_times.values())}"
        s += ")"
        return s

    def __hash__(self):
        tmp = (self.t, self.req, self.req_times, self.D, self.E)
        return hash(tuple(tuple(x) for x in tmp))


@dataclass(frozen=True, repr=False)
class Action:
    sat_id: int  # Satellite id
    ts: float  # Start time
    dur: float  # Duration
    req: Request  # Request

    def __repr__(self):
        a = "Action("
        a += f"sat={self.sat_id}, "
        a += f"usr={self.req.usr_id if self.req else None}, "
        a += f"req={self.req.id if self.req else None}, "
        a += f"ts={self.ts:.2f}, "
        a += f"dur={self.dur:.2f}"
        a += ")"
        return a


def logging_decorator(func):
    def wrapper(*args, **kwargs):
        logging.debug(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper


# *************************************************************************************************
# Problem Data
# *************************************************************************************************


@dataclass
class PntSchedulingProblem:
    # Scenario
    t_step: Time
    t_final: Time
    requests: list[Request]
    service_windows: list[ServiceWindow]
    transition_times: np.ndarray  # (N_req x N_req)
    CN0: np.ndarray  # (N_sat x N_users x N_steps)

    # Resources
    min_energy: Energy
    max_energy: Energy
    min_data: Data
    max_data: Data
    payload_energy_gen: float
    payload_data_gen: float
    energy_gen: np.ndarray  # (N_sat x N_steps)
    data_gen: np.ndarray  # (N_sat x N_steps)

    # Online
    current_policy: list[tuple[State, Action]] = None
    current_time: float = 0
    constr: bool = False
    constr_windows: list[ServiceWindow] = None
    constr_windows_dict: dict[ReqId, list[ServiceWindow]] = None
    required_resources: dict[SatId, tuple[Time, Data, Energy]] = None
    action_starts: np.ndarray = None
    action_end: np.ndarray = None
    action_sats_ids: np.ndarray = None

    # Computed
    service_windows_dict: dict[ReqId, list[ServiceWindow]] = None
    CN0_norm: np.ndarray = None  # Normalized CN0
    N_sat: int = None  # Number of satellites
    times: np.ndarray = None  # Time vector

    def __post_init__(self):
        self.service_windows_dict = {
            req.id: [win for win in self.service_windows if win.usr_id == req.usr_id]
            for req in self.requests
        }
        self.CN0_norm = self.CN0 / np.nanmax(self.CN0, axis=2)[:, :, None]
        self.N_sat = len(set(w.sat_id for w in self.service_windows))
        self.times = np.arange(0, self.t_final + self.t_step, self.t_step)

    # *********************************************************************************************
    # SMDP Functions
    # *********************************************************************************************

    def available_actions(self, s: State, N_max: int, dur_min: float) -> list[Action]:
        """
        Get available actions for a given state

        Args:
            s (State): Current state
            N_max (int): Maximum number of actions
            T_min (float): Minimum action duration

        Returns:
            list[Action]: List of available actions
        """

        # Satellite to select an action (min ti, and min id)
        sat_id = np.where(s.t == np.min(s.t))[0][0]
        if s.t[sat_id] >= self.t_final:
            return []

        ti = s.t[sat_id]
        Ei = s.E[sat_id]
        Di = s.D[sat_id]
        req_id = s.req[sat_id].id if s.req[sat_id] else None

        requests = [
            req
            for req in self.requests
            if req.ta <= self.current_time  # Arrived
            and ti < req.te  # Still valid
            and s.req_times[req.id] < req.dur  # Not fulfilled
        ]

        actions = [
            Action(
                sat_id=sat_id,
                ts=s.t[sat_id],
                dur=min(dur_min, self.t_final - s.t[sat_id]),
                req=None,
            )
        ]

        for req in requests:
            dur = min(dur_min, req.dur - s.req_times[req.id])

            # Transition time
            t_trans = self.transition_times[sat_id, req_id, req.id] if req_id else ti

            # User available
            sat_id_other = [
                sat_id_
                for sat_id_ in range(self.N_sat)
                if s.req[sat_id_] and s.req[sat_id_].usr_id == req.usr_id
            ]
            t_user = s.t[sat_id_other[0]] if sat_id_other else ti

            # Service start
            t_earliest = max(t_trans, t_user)

            # Windows
            win_dict = (
                self.service_windows_dict[req.id]
                if not (self.constr and self.current_policy)
                else self.constr_windows_dict[req.id]
            )
            win_list = [
                win
                for win in win_dict
                if win.sat_id == sat_id  # Right satellite
                and win.te >= t_earliest + min(dur, win.te - win.ts)
                # Clipped for short windows below T_min
            ]

            # Energy and data
            actions_req = []
            for win in win_list:
                dur_win = min(dur, win.te - win.ts)
                t_win = max(t_earliest, win.ts)
                # Earliest of action
                Es = Ei + self.energy_gen_func(sat_id, ti, t_win)
                Ds = Di + self.data_gen_func(sat_id, ti, t_win)
                # Required energy and data
                Ea = self.payload_energy_gen * dur_win
                Da = self.payload_data_gen * dur_win

                t_E = self.time_for_energy(sat_id, t_win, self.min_energy - Ea - Es)
                t_D = self.time_for_data(sat_id, t_win, self.max_data - Da - Ds)

                if t_E is None or t_D is None:
                    # Not enough resources
                    continue

                # Subtract the time of the action
                t_E -= dur_win
                t_D -= dur_win

                ts = max(t_win, t_E, t_D, win.ts)
                Es = Ei + self.energy_gen_func(sat_id, ti, ts)
                Ds = Di + self.data_gen_func(sat_id, ti, ts)

                # Check if the action is feasible
                if ts + dur_win > win.te:
                    continue

                # Check required resources after the action
                if self.constr and self.current_policy:
                    constraints = [
                        const
                        for const in self.required_resources[sat_id]
                        if const[0] >= ts + dur_win
                    ]
                    if constraints:
                        tc, Ec, Dc = constraints[0]
                        if (
                            Es + Ea + self.energy_gen_func(sat_id, ts, tc) < Ec
                            or Ds + Da + self.data_gen_func(sat_id, ts, tc) > Dc
                        ):
                            continue

                a = Action(sat_id=sat_id, ts=ts, dur=dur_win, req=req)
                self.transition_function(s, a)
                actions_req.append(a)

            actions.extend(actions_req)

        actions.sort(key=lambda a: a.ts)
        actions = actions[:N_max]
        return actions

    def initial_state(self) -> State:
        """
        Get initial state

        Returns:
            State: Initial state
        """
        t = self.get_current_policy_index()

        if t is None:
            return State(
                t=[0.0 for _ in range(self.N_sat)],
                req=[None for _ in range(self.N_sat)],
                req_times={req.id: 0.0 for req in self.requests},
                D=[(self.min_data + self.max_data) / 2.0 for _ in range(self.N_sat)],
                E=[
                    (self.min_energy + self.max_energy) / 2.0 for _ in range(self.N_sat)
                ],
            )

        s, _ = deepcopy(self.current_policy[t])
        for sat_id in range(self.N_sat):
            if s.t[sat_id] < self.current_time:
                s.t[sat_id] = self.current_time

            if self.constr and self.current_policy:
                # Account for the current schedule
                sf = self.current_policy[-1][0]
                for req_id in s.req_times:
                    s.req_times[req_id] = sf.req_times[req_id]
        return s

    def transition_function(self, s: State, a: Action) -> State:
        """
        Transition function

        Args:
            s (State): Current state
            a (Action): Action

        Returns:
            State: New state
        """

        time = a.ts + a.dur
        sat_id = a.sat_id

        if a.req:
            usr_id = a.req.usr_id
            req_id = a.req.id
            payload_on = usr_id >= 0

            assert a.dur <= a.req.dur - s.req_times[req_id]
            duration = min(a.dur, a.req.dur - s.req_times[req_id])

        else:
            payload_on = False
            duration = a.dur

        data = max(
            self.min_data,
            s.D[sat_id]
            + self.payload_data_gen * duration * payload_on
            + self.data_gen_func(sat_id, s.t[sat_id], time),
        )
        energy = min(
            self.max_energy,
            s.E[sat_id]
            + self.payload_energy_gen * duration * payload_on
            + self.energy_gen_func(sat_id, s.t[sat_id], time),
        )

        if data > self.max_data + ABS_TOL:
            tmp = max(
                self.min_data,
                s.D[sat_id]
                + self.payload_data_gen * duration * payload_on
                + self.data_gen_func(sat_id, s.t[sat_id], time),
            )

        assert s.t[sat_id] <= a.ts + ABS_TOL
        assert data >= self.min_data - ABS_TOL
        assert data <= self.max_data + ABS_TOL
        assert energy <= self.max_energy + ABS_TOL
        assert energy >= self.min_energy - ABS_TOL

        sp = deepcopy(s)
        sp.t[sat_id] = time
        sp.req[sat_id] = a.req
        if a.req:
            sp.req_times[req_id] = s.req_times[req_id] + duration
            assert sp.req_times[req_id] <= a.req.dur
        sp.D[sat_id] = data
        sp.E[sat_id] = energy
        return sp

    def integrate_normalized_CN0(self, ts, te, request_id, sat_id) -> float:
        """
        Integrate normalized CN0

        Args:
            ts (float): Start time
            te (float): End time
            request_id (int): Request id
            sat_id (int): Satellite id

        Returns:
            float: Integrated normalized CN0
        """

        i_s, i_e = self.get_discrete_index(ts, te)
        cn0 = self.CN0_norm[sat_id, request_id, i_s:i_e]
        cn0[np.isnan(cn0)] = 0
        return np.sum(cn0) * self.t_step

    def reward_function(self, s: State, a: Action) -> float:
        """
        Reward function

        Args:
            s (State): Current state
            a (Action): Action

        Returns:
            float: Reward
        """
        if a.req is None:
            return 0

        sat_id = a.sat_id
        same_req = int(a.req and s.req[sat_id] == a.req and s.t[sat_id] == a.ts)

        payload_on = a.req.usr_id >= 0
        if payload_on:
            cn0 = self.integrate_normalized_CN0(
                a.ts, a.ts + a.dur, a.req.usr_id, a.sat_id
            )
            bonus = same_req * 1.0
            mult = 0.9 ** (a.ts - a.req.ts)
            return mult * (cn0 + bonus)
        else:
            return 0

    # *********************************************************************************************
    # Resources
    # *********************************************************************************************

    def get_discrete_index(self, ts, te):
        if isinstance(ts, np.ndarray):
            i_s = (ts / self.t_step).astype(int)
        else:
            i_s = int(ts / self.t_step)

        if isinstance(te, np.ndarray):
            i_e = (te / self.t_step).astype(int)
        else:
            i_e = int(te / self.t_step)

        return i_s, i_e

    def intersection_with_current_actions(self, sat_id, ts, te):
        if not (self.constr and self.current_policy):
            return 0
        starts = np.maximum(ts, self.action_starts)
        ends = np.minimum(te, self.action_ends)
        same_sat = (self.action_sats_ids == sat_id).astype(int)
        return np.sum(np.maximum(0, (ends - starts) * same_sat))

    def data_gen_func(self, sat_id, ts, te, constr=True):
        i_s, i_e = self.get_discrete_index(ts, te)
        duration = self.intersection_with_current_actions(sat_id, ts, te) * int(constr)
        return np.sum(self.data_gen[sat_id, i_s:i_e]) + self.payload_data_gen * duration

    def energy_gen_func(self, sat_id, ts, te, constr=True):
        i_s, i_e = self.get_discrete_index(ts, te)
        duration = self.intersection_with_current_actions(sat_id, ts, te) * int(constr)
        return (
            np.sum(self.energy_gen[sat_id, i_s:i_e])
            + self.payload_energy_gen * duration
        )

    def time_for_energy(self, sat_id, ts, energy):
        i_s = int(ts / self.t_step)
        for i_e in range(i_s, len(self.energy_gen[sat_id])):
            te = i_e * self.t_step
            if self.energy_gen_func(sat_id, ts, te) >= energy - ABS_TOL:
                return te

    def time_for_data(self, sat_id, ts, data):
        i_s = int(ts / self.t_step)
        for i_e in range(i_s, len(self.data_gen[sat_id])):
            te = i_e * self.t_step
            if self.data_gen_func(sat_id, ts, te) <= data + ABS_TOL:
                return te

    # *********************************************************************************************
    # Constraints
    # *********************************************************************************************

    def get_current_policy_index(self):
        if self.current_policy is None:
            return None
        for t, (s, a) in enumerate(self.current_policy):
            if (
                np.min(s.t) > self.current_time  # Minimum time for all satellites
                or a is None  # No more actions
                or a.ts > self.current_time  # Next action is past the current time
                or (a.req is None and a.ts + a.dur > self.current_time)
            ):
                return t
        return len(self.current_policy) - 1

    def set_current_time(self, current_time):
        self.current_time = current_time

    def set_contrained(self, constr):
        self.constr = constr

    def set_current_policy(self, policy: list[tuple[State, Action]]):
        self.current_policy = deepcopy(policy)

        if not (self.constr and self.current_policy):
            self.constr_windows = None
            self.constr_windows_dict = None
            self.action_starts = None
            self.action_ends = None
            self.action_sats_ids = None
            self.required_resources = None
            return

        # Create new windows for the constr problem
        actions_by_win = {
            win.id: [
                a
                for _, a in self.current_policy
                if a and a.req and (a.ts + a.dur > win.ts and a.ts < win.te)
            ]
            for win in self.service_windows
        }

        win_list = deepcopy(self.service_windows)
        actions = [a for _, a in self.current_policy if a]
        for a in actions:
            new_win_list = []
            for win in win_list:
                if a.ts + a.dur <= win.ts or a.ts >= win.te:
                    # No intersection
                    new_win_list.append(win)

                elif a.ts <= win.ts and a.ts + a.dur >= win.te:
                    # Fully covered
                    pass

                elif a.ts > win.ts and a.ts + a.dur < win.te:
                    # Middle part
                    win_ = ServiceWindow(win.sat_id, win.usr_id, win.ts, a.ts)
                    win_list.append(win_)
                    win_ = ServiceWindow(win.sat_id, win.usr_id, a.ts + a.dur, win.te)
                    win_list.append(win_)

                elif a.ts <= win.ts:
                    # Left part
                    win_ = ServiceWindow(win.sat_id, win.usr_id, a.ts + a.dur, win.te)
                    win_list.append(win_)

                elif a.ts + a.dur >= win.te:
                    # Right part
                    win_ = ServiceWindow(win.sat_id, win.usr_id, win.ts, a.ts)
                    win_list.append(win_)
                else:
                    raise ValueError("Invalid case")

            win_list = new_win_list

        self.constr_windows = win_list
        self.constr_windows_dict = {
            req.id: [win for win in win_list if win.usr_id == req.usr_id]
            for req in self.requests
        }
        self.action_starts = np.array([a.ts for _, a in self.current_policy if a])
        self.action_ends = np.array([a.ts + a.dur for _, a in self.current_policy if a])
        self.action_sats_ids = np.array([a.sat_id for _, a in self.current_policy if a])

        self.compute_required_resources(policy)

    def compute_required_resources(self, policy: list[tuple[State, Action]]):
        required_resources = dict[SatId, list[tuple[Time, Energy, Data]]]()
        for sat_id in range(self.N_sat):
            required_resources[sat_id] = list[tuple[Time, Energy, Data]]()
            actions = [a for s, a in policy[::-1] if a and a.req and a.sat_id == sat_id]
            if len(actions) == 0:
                continue
            d = self.max_data
            e = self.min_energy
            a = actions[0]
            t_prev = a.ts + a.dur
            for a in actions:
                # From end of current action to start of next action
                d -= self.data_gen_func(sat_id, a.ts + a.dur, t_prev)
                e -= self.energy_gen_func(sat_id, a.ts + a.dur, t_prev)

                # Clip values
                d = min(d, self.max_data)
                e = max(e, self.min_energy)

                # From start to end of current action
                d -= self.data_gen_func(sat_id, a.ts, a.ts + a.dur)
                e -= self.energy_gen_func(sat_id, a.ts, a.ts + a.dur)
                d -= self.payload_data_gen * a.dur
                e -= self.payload_energy_gen * a.dur

                required_resources[sat_id].append((a.ts, d, e))
                t_prev = a.ts

            # Sort by time
            required_resources[sat_id].sort(key=lambda x: x[0])

        self.required_resources = required_resources

    def merge_policies(
        self,
        policy: list[tuple[State, Action]],
    ) -> list[tuple[State, Action]]:
        current_time = self.current_time
        constr = self.constr
        index = (
            -1
            if (self.constr and self.current_policy)
            else self.get_current_policy_index()
        )

        self.current_time = 0
        self.constr = False

        actions = [
            a for _, a in self.current_policy[:index] if a and self.current_policy
        ]
        actions.extend([a for _, a in policy if a])
        actions.sort(key=lambda a: a.ts)
        new_policy = []
        s = self.current_policy[0][0] if self.current_policy else policy[0][0]
        for a in actions:
            if a.req:
                new_policy.append((s, a))
                s = self.transition_function(s, a)
        new_policy.append((s, None))

        self.current_time = current_time
        self.constr = constr

        return new_policy

    # *********************************************************************************************
    # Utils
    # *********************************************************************************************

    def get_arrival_times(self):
        return sorted(list(set(req.ta for req in self.requests)))

    def total_reward(self, policy: list[tuple[State, Action]], gamma: float) -> float:
        """
        Calculate the total reward of a policy

        Args:
            policy (list[tuple[State, Action]]): Policy
            gamma (float): Discount factor

        Returns:
            float: Total reward
        """

        rewards = np.array([self.reward_function(s, a) for s, a in policy[:-1]])
        discounts = np.array([gamma ** (a.ts + 0.5 * a.dur) for _, a in policy[:-1]])
        return np.sum(rewards * discounts)

    def percentage_completed(self, policy: list[tuple[State, Action]]) -> list[float]:
        """
        Calculate the percentage of completed requests

        Args:
            policy (list[tuple[State, Action]]): Policy

        Returns:
            list[float]: Percentage of completed requests
        """

        s0 = policy[0][0]
        sf = policy[-1][0]
        return [
            (
                sf.req_times[req.id] / req.dur * 100
                if req.ta <= self.current_time
                else np.NaN
            )
            for req in self.requests
            if req.id >= 0
        ]

    def duration_fullfilled(self, policy: list[tuple[State, Action]]) -> float:
        """
        Calculate the percentage of completed requests

        Args:
            policy (list[tuple[State, Action]]): Policy

        Returns:
            float: Percentage of completed requests
        """

        s0 = policy[0][0]
        sf = policy[-1][0]
        return (
            sum(sf.req_times[req.id] for req in self.requests if req.id >= 0)
            / sum(req.dur for req in self.requests if req.id >= 0)
            * 100
        )

    def clean_policy(
        self, policy: list[tuple[State, Action]]
    ) -> list[tuple[State, Action]]:
        """
        Clean policy by removing consecutive actions

        Args:
            policy (list[tuple[State, Action]]): Policy

        Returns:
            list[tuple[State, Action]]: Cleaned policy
        """

        actions_by_sat = [list[Action]() for _ in range(self.N_sat)]
        for _, a in policy[:-1]:
            sat_id = a.sat_id
            last_a = actions_by_sat[sat_id][-1] if actions_by_sat[sat_id] else None
            if (
                not last_a  # Empty list
                or a.req != last_a.req  # Different request
                or a.ts != last_a.ts + last_a.dur  # Not consecutive
            ):
                actions_by_sat[sat_id].append(a)
            else:
                new_a = Action(
                    sat_id=a.sat_id, ts=last_a.ts, dur=last_a.dur + a.dur, req=a.req
                )
                actions_by_sat[sat_id][-1] = new_a
        all_actions = sorted(
            [a for actions in actions_by_sat for a in actions], key=lambda a: a.ts
        )

        new_policy = []
        s = policy[0][0]
        for a in all_actions:
            if a.req:
                new_policy.append((s, a))
                s = self.transition_function(s, a)
        new_policy.append((s, None))
        return new_policy
