import numpy as np
import logging
from dataclasses import dataclass, field
from copy import deepcopy
from typing import ClassVar

from dataclasses import dataclass


def reset_id_counters():
    Request._next_id = 0
    ServiceWindow._next_id = 0


class UsrId(int):
    pass


class SatId(int):
    pass


class ReqId(int):
    pass


class Data(float):
    pass


class Energy(float):
    pass


class Time(float):
    pass


@dataclass(frozen=True, repr=True)
class Request:
    user_id: UsrId  # User id
    start: float  # Start time
    end: float  # End time
    duration: float  # Duration
    priority: int = 1  # Priority
    arrival: float = 0  # Arrival time
    id: int = field(default_factory=lambda: Request._next_id)  # Request id

    _next_id: ClassVar[int] = 0

    def __post_init__(self):
        Request._next_id += 1


@dataclass(frozen=True, repr=True)
class ServiceWindow:
    satellite_id: SatId  # Satellite id
    user_id: int  # User id
    start: float  # Start time
    end: float  # End time
    id: int = field(default_factory=lambda: ServiceWindow._next_id)  # Window id

    _next_id: ClassVar[int] = 0

    def __post_init__(self):
        ServiceWindow._next_id += 1


@dataclass(frozen=True, repr=False)
class State:
    times: list[float]  # Current time
    requests: list[Request]  # Last request
    request_times: dict[ReqId, float]  # Request time
    data: list[float]  # Data onboard
    energy: list[float]  # Energy onboard

    def __repr__(self):
        s = "State("
        s += f"t={np.round(self.times, 2)}, "
        s += f"d={np.round(self.data, 2)}, "
        s += f"e={np.round(self.energy, 2)}, "
        s += f"req={list(self.request_times.values())}"
        s += ")"
        return s

    def __hash__(self):
        return hash(
            (
                tuple(self.times),
                tuple(self.requests),
                tuple(self.request_times.items()),
                tuple(self.data),
                tuple(self.energy),
            )
        )


@dataclass(frozen=True, repr=False)
class Action:
    satellite_id: int  # Satellite id
    start: float  # Start time
    duration: float  # Duration
    request: Request  # Request

    def __repr__(self):
        a = "Action("
        a += f"sat={self.satellite_id}, "
        a += f"usr={self.request.user_id if self.request is not None else None}, "
        a += f"req={self.request.id if self.request is not None else None}, "
        a += f"start={self.start:.2f}, "
        a += f"dur={self.duration:.2f}"
        a += ")"
        return a


def logging_decorator(func):
    def wrapper(*args, **kwargs):
        logging.debug(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper


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
    constrained: bool = False
    constrained_windows: list[ServiceWindow] = None
    constrained_windows_dict: dict[ReqId, list[ServiceWindow]] = None
    required_resources: dict[SatId, tuple[Time, Data, Energy]] = None

    # Computed
    service_windows_dict: dict[ReqId, list[ServiceWindow]] = None
    CN0_norm: np.ndarray = None  # Normalized CN0
    N_sat: int = None  # Number of satellites
    times: np.ndarray = None  # Time vector

    def __post_init__(self):
        self.service_windows_dict = {
            req.id: [win for win in self.service_windows if win.user_id == req.user_id]
            for req in self.requests
        }
        self.CN0_norm = self.CN0 / np.nanmax(self.CN0, axis=2)[:, :, None]
        self.N_sat = len(set(w.satellite_id for w in self.service_windows))
        self.times = np.arange(0, self.t_final + self.t_step, self.t_step)

    def set_current_time(self, current_time):
        self.current_time = current_time

    def set_current_policy(
        self, policy: list[tuple[State, Action]], constrained: bool = False
    ):
        self.current_policy = deepcopy(policy)
        self.constrained = constrained

        if not constrained:
            return

        # Create new windows for the constrained problem
        actions_by_win = {
            win.id: [
                a
                for _, a in self.current_policy
                # Not last action
                if a is not None
                # Correct satellite
                and a.satellite_id == win.satellite_id
                # Correct window
                and (a.start >= win.start and a.start + a.duration <= win.end)
            ]
            for win in self.service_windows
        }

        win_list = list[ServiceWindow]()
        for win_id, actions in actions_by_win.items():
            win = self.service_windows[win_id]

            ts = win.start
            for a in actions:
                te = a.start
                if te > ts:
                    win_list.append(
                        ServiceWindow(
                            satellite_id=win.satellite_id,
                            user_id=win.user_id,
                            start=ts,
                            end=te,
                        )
                    )
                ts = a.start + a.duration

            if ts < win.end:
                win_list.append(
                    ServiceWindow(
                        satellite_id=win.satellite_id,
                        user_id=win.user_id,
                        start=ts,
                        end=win.end,
                    )
                )

        self.constrained_windows = win_list
        self.constrained_windows_dict = {
            req.id: [win for win in win_list if win.user_id == req.user_id]
            for req in self.requests
        }

    def data_gen_func(self, sat_id, ts, te):
        i_s, i_e = self.get_discrete_index(ts, te)
        if isinstance(i_s, np.ndarray) or isinstance(i_e, np.ndarray):
            return np.array(
                [
                    np.sum(self.data_gen[sat_id, i_s_:i_e_])
                    for i_s_, i_e_ in zip(*np.broadcast_arrays(i_s, i_e))
                ]
            )
        else:
            return np.sum(self.data_gen[sat_id, i_s:i_e])

    def energy_gen_func(self, sat_id, ts, te):
        i_s, i_e = self.get_discrete_index(ts, te)
        if isinstance(i_s, np.ndarray) or isinstance(i_e, np.ndarray):
            return np.array(
                [
                    np.sum(self.energy_gen[sat_id, i_s_:i_e_])
                    for i_s_, i_e_ in zip(*np.broadcast_arrays(i_s, i_e))
                ]
            )
        else:
            return np.sum(self.energy_gen[sat_id, i_s:i_e])

    def time_for_energy(self, sat_id, ts, energy):
        i_s = int(ts / self.t_step)
        di = np.where(np.cumsum(self.energy_gen[sat_id, i_s:]) >= energy)[0]
        if len(di) == 0:
            return None
        return self.times[i_s + di[0]]

    def time_for_data(self, sat_id, ts, data):
        i_s = int(ts / self.t_step)
        di = np.where(np.cumsum(self.data_gen[sat_id, i_s:]) <= data)[0]
        if len(di) == 0:
            return None
        return self.times[i_s + di[0]]

    def get_arrival_times(self):
        return sorted(list(set(req.arrival for req in self.requests)))

    def compute_required_resources(self, policy: list[tuple[State, Action]]):
        required_resources = dict[SatId, list[tuple[Time, Energy, Data]]]()
        for sat_id in range(self.N_sat):
            required_resources[sat_id] = list[tuple[Time, Energy, Data]]()
            actions = [
                a
                for s, a in policy[::-1]
                if a is not None and a.request is not None and a.satellite_id == sat_id
            ]
            if len(actions) == 0:
                continue
            d = self.max_data
            e = self.min_energy
            a = actions[0]
            t_prev = a.start + a.duration
            for a in actions:

                # From end of current action to start of next action
                d -= self.data_gen_func(a.start + a.duration, t_prev)
                e -= self.energy_gen_func(a.start + a.duration, t_prev)

                # Clip values
                d = min(d, self.max_data)
                e = max(e, self.min_energy)

                # From start to end of current action
                d -= self.data_gen_func(a.start, a.start + a.duration)
                e -= self.energy_gen_func(a.start, a.start + a.duration)
                d -= self.payload_data_gen * a.duration
                e -= self.payload_energy_gen * a.duration

                required_resources[sat_id].append((a.start, d, e))
                t_prev = a.start

            # Sort by time
            required_resources[sat_id].sort(key=lambda x: x[0])

        self.required_resources = required_resources

    def available_actions(self, s: State, N_max: int, T_min: float) -> list[Action]:
        """
        Get available actions for a given state

        Args:
            s (State): Current state
            N_max (int): Maximum number of actions
            T_min (float): Minimum action duration

        Returns:
            list[Action]: List of available actions
        """

        # Satellite to select an action
        sat_id = np.argmin(s.times)
        if s.times[sat_id] >= self.t_final:
            return []

        ti = s.times[sat_id]
        Ei = s.energy[sat_id]
        Di = s.data[sat_id]
        req_id = s.requests[sat_id].id if s.requests[sat_id] is not None else None

        requests = [
            req
            for req in self.requests
            if req.arrival <= self.current_time  # Arrived
            and ti < req.end  # Still valid
            and s.request_times[req.id] < req.duration  # Not fulfilled
        ]

        actions = [
            Action(
                satellite_id=sat_id,
                start=s.times[sat_id],
                duration=min(T_min, self.t_final - s.times[sat_id]),
                request=None,
            )
        ]
        for req in requests:
            duration = min(T_min, req.duration - s.request_times[req.id])

            # Transition time
            t_trans = (
                self.transition_times[sat_id, req_id, req.id]
                if req_id is not None
                else ti
            )

            # User available
            sat_id_other = [
                sat_id_
                for sat_id_ in range(self.N_sat)
                if s.requests[sat_id_] is not None
                and s.requests[sat_id_].user_id == req.user_id
            ]
            t_user = s.times[sat_id_other[0]] if sat_id_other else ti

            # Service start
            t_earliest = max(t_trans, t_user)

            # Windows
            win_list = [
                win
                for win in self.service_windows_dict[req.id]
                if win.satellite_id == sat_id  # Right satellite
                and win.end
                >= t_earliest + min(duration, win.end - win.start)  # Enough time
                # Clipped for short windows below T_min
            ]
            t_win_list = [max(t_earliest, win.start) for win in win_list]

            # Energy and data
            actions_req = []
            for win in win_list:
                T = min(duration, win.end - win.start)
                E = self.payload_energy_gen * T
                D = self.payload_data_gen * T
                t_E = self.time_for_energy(sat_id, t_earliest, self.min_energy - Ei - E)
                t_D = self.time_for_data(sat_id, t_earliest, self.max_data - Di - D)
                ts = max(t_earliest, t_E, t_D, win.start)
                if ts + T <= win.end:
                    actions_req.append(
                        Action(satellite_id=sat_id, start=ts, duration=T, request=req)
                    )

            actions.extend(actions_req)

        return actions

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

    def transition_function(self, s: State, a: Action) -> State:
        """
        Transition function

        Args:
            s (State): Current state
            a (Action): Action

        Returns:
            State: New state
        """

        time = a.start + a.duration
        sat_id = a.satellite_id

        if a.request is not None:
            usr_id = a.request.user_id
            req_id = a.request.id
            payload_on = usr_id >= 0

            assert a.duration <= a.request.duration - s.request_times[req_id]
            duration = min(
                a.duration,
                a.request.duration - s.request_times[req_id],
            )

        else:
            payload_on = False
            duration = a.duration

        data = max(
            self.min_data,
            s.data[sat_id]
            + self.payload_data_gen * duration * payload_on
            + self.data_gen_func(sat_id, s.times[sat_id], time),
        )
        energy = min(
            self.max_energy,
            s.energy[sat_id]
            + self.payload_energy_gen * duration * payload_on
            + self.energy_gen_func(sat_id, s.times[sat_id], time),
        )

        eps = 1e-6
        assert s.times[sat_id] <= a.start + eps
        assert data >= self.min_data - eps
        assert data <= self.max_data + eps
        assert energy <= self.max_energy + eps
        assert energy >= self.min_energy - eps

        sp = deepcopy(s)
        sp.times[sat_id] = time
        sp.requests[sat_id] = a.request
        if a.request is not None:
            sp.request_times[req_id] = s.request_times[req_id] + duration
            assert sp.request_times[req_id] <= a.request.duration
        sp.data[sat_id] = data
        sp.energy[sat_id] = energy
        return sp

    def integrate_normalized_CN0(self, ts, te, request_id, satellite_id) -> float:
        """
        Integrate normalized CN0

        Args:
            ts (float): Start time
            te (float): End time
            request_id (int): Request id
            satellite_id (int): Satellite id

        Returns:
            float: Integrated normalized CN0
        """

        i_s, i_e = self.get_discrete_index(ts, te)
        cn0 = self.CN0_norm[satellite_id, request_id, i_s:i_e]
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
        if a.request is None:
            return 0

        sat_id = a.satellite_id
        same_req = int(
            a.request is not None
            and s.requests[sat_id] == a.request
            and s.times[sat_id] == a.start
        )

        payload_on = a.request.user_id >= 0
        if payload_on:
            cn0 = self.integrate_normalized_CN0(
                a.start, a.start + a.duration, a.request.user_id, a.satellite_id
            )
            bonus = same_req * 1.0
            mult = 0.9 ** (a.start - a.request.start)
            return mult * (cn0 + bonus)
        else:
            return 0

    def get_current_policy_index(self):
        if self.current_policy is None:
            return None
        for t, (s, a) in enumerate(self.current_policy):
            if (
                # All satellites have received requests
                np.min(s.times) > self.current_time
                # No more actions
                or a is None
                # Next action is past the current time
                or a.start > self.current_time
                # Next action is not a request
                or (a.request is None and a.start + a.duration > self.current_time)
            ):
                return t
        return len(self.current_policy) - 1

    def initial_state(self) -> State:
        """
        Get initial state

        Returns:
            State: Initial state
        """
        t = self.get_current_policy_index()

        if t is None:
            return State(
                times=[0.0 for _ in range(self.N_sat)],
                requests=[None for _ in range(self.N_sat)],
                request_times={req.id: 0.0 for req in self.requests},
                data=[(self.min_data + self.max_data) / 2.0 for _ in range(self.N_sat)],
                energy=[
                    (self.min_energy + self.max_energy) / 2.0 for _ in range(self.N_sat)
                ],
            )

        s, _ = self.current_policy[t]
        for sat_id in range(self.N_sat):
            if s.times[sat_id] < self.current_time:
                s.times[sat_id] = self.current_time

            if self.constrained:
                # Account for the current schedule
                sf = self.current_policy[-1][0]
                for req_id in s.request_times:
                    s.request_times[req_id] = sf.request_times[req_id]
        return s

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
        discounts = np.array(
            [gamma ** (a.start + 0.5 * a.duration) for _, a in policy[:-1]]
        )
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
                sf.request_times[req.id] / req.duration * 100
                if req.arrival <= self.current_time
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
            sum(sf.request_times[req.id] for req in self.requests if req.id >= 0)
            / sum(req.duration for req in self.requests if req.id >= 0)
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
            sat_id = a.satellite_id
            last_a = actions_by_sat[sat_id][-1] if actions_by_sat[sat_id] else None
            if (
                not last_a  # Empty list
                or a.request != last_a.request  # Different request
                or a.start != last_a.start + last_a.duration  # Not consecutive
            ):
                actions_by_sat[sat_id].append(a)
            else:
                new_a = Action(
                    satellite_id=a.satellite_id,
                    start=last_a.start,
                    duration=last_a.duration + a.duration,
                    request=a.request,
                )
                actions_by_sat[sat_id][-1] = new_a
        all_actions = sorted(
            [a for actions in actions_by_sat for a in actions], key=lambda a: a.start
        )

        new_policy = []
        s = policy[0][0]
        for a in all_actions:
            if a.request is not None:
                new_policy.append((s, a))
                s = self.transition_function(s, a)
        new_policy.append((s, None))
        return new_policy
