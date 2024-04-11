import numpy as np
import logging
from dataclasses import dataclass, field
from copy import deepcopy
from typing import ClassVar


def reset_id_counters():
    Request._next_id = 0
    ServiceWindow._next_id = 0


@dataclass(frozen=True, repr=True)
class Request:
    user_id: int  # User id
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
    satellite_id: int  # Satellite id
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
    request_times: dict[int, float]  # Request time
    data: list[float]  # Data onboard
    energy: list[float]  # Energy onboard

    def __repr__(self):
        s = "State("
        s += f"t={np.round(self.times, 2)}, "
        s += f"d={np.round(self.data, 2)}, "
        s += f"e={np.round(self.energy, 2)}, "
        s += f"request_time={self.request_times}"
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
        a += f"user={self.request.user_id if self.request is not None else None}, "
        a += f"request={self.request.id if self.request is not None else None}, "
        a += f"start={self.start:.2f}, "
        a += f"duration={self.duration:.2f}"
        a += ")"
        return a


def logging_decorator(func):
    def wrapper(*args, **kwargs):
        logging.debug(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper


@dataclass()
class PntSchedulingProblem:

    def __init__(
        self,
        # Problem data
        time_step: float,  # Time step
        requests: list[Request],  # List of requests
        service_windows: list[ServiceWindow],  # List of service windows
        transition_times: np.ndarray,  # Transition times matrix
        CN0,  # Carrier-to-noise density ratio [dB-Hz]
        # Satelite
        min_energy: float,  # Minimum energy
        max_energy: float,  # Maximum energy
        min_data: float,  # Minimum data
        max_data: float,  # Maximum data
        payload_power_gen: float,  # Payload power generation
        payload_data_gen: float,  # Payload data generation
        energy_gen_func: callable,  # Energy generation function
        data_gen_func: callable,  # Data generation function
    ):
        self.time_step = time_step
        self.requests = sorted(requests, key=lambda req: req.start)
        self.service_windows = service_windows
        self.service_windows_dict = {
            req.id: [win for win in service_windows if win.user_id == req.user_id]
            for req in requests
        }
        self.transition_times = transition_times
        self.CN0_norm = CN0 / np.nanmax(CN0, axis=2)[:, :, None]
        self.N_satellites = len(set(w.satellite_id for w in service_windows))

        self.min_energy = min_energy
        self.max_energy = max_energy
        self.min_data = min_data
        self.max_data = max_data
        self.payload_power_gen = payload_power_gen
        self.payload_data_gen = payload_data_gen
        self.energy_gen_func = energy_gen_func
        self.data_gen_func = data_gen_func
        self.tf = max([r.end for r in requests])

    def available_actions(self, s: State, N_max: int, d_min: float) -> list[Action]:
        """
        Get available actions for a given state

        Args:
            s (State): Current state
            N_max (int): Maximum number of actions
            d_min (float): Minimum action duration

        Returns:
            list[Action]: List of available actions
        """

        # Satellite to select an action
        sat_id = np.argmin(s.times)
        if s.times[sat_id] >= self.tf:
            return []

        # List of users being served by other satellites
        other_sat_user_ids = [
            req.user_id
            for other_sat_id, req in enumerate(s.requests)
            if req is not None and other_sat_id != sat_id
        ]

        # List of available windows for each request
        available_windows = {
            req.id: [
                win
                for win in self.service_windows_dict[req.id]
                if win.satellite_id == sat_id and win.end >= s.times[sat_id]
            ]
            for req in self.requests
        }

        # List of available requests
        available_requests = [
            req
            for req in self.requests
            if
            # There is still time to fulfill the request
            req.end >= s.times[sat_id]
            # The request is not completed
            and s.request_times[req.id] < req.duration
            # There is a service window for the request
            and len(available_windows[req.id]) > 0
            # The user is not being served by another satellite
            and req.user_id not in other_sat_user_ids
        ]

        actions = [
            Action(
                satellite_id=sat_id,
                start=s.times[sat_id],
                duration=d_min,
                request=None,
            )
        ]
        actions_count = 0
        for req in available_requests:
            # Transition time
            trans_time = (
                0
                if req.user_id < 0 or s.requests[sat_id] is None
                else self.transition_times[
                    sat_id, s.requests[sat_id].user_id, req.user_id
                ]
            )

            for win in available_windows[req.id]:  # Only consider the first window
                # Action start and duration
                a_start = max(
                    win.start,
                    s.times[sat_id] + trans_time,
                )
                a_duration = min(
                    d_min,  # Minimum duration (parameter)
                    win.end - a_start,  # Window end
                    req.duration - s.request_times[req.user_id],  # Requested duration
                )
                if a_duration <= 0:
                    continue

                # Resources
                energy_gen = self.energy_gen_func(s.times[sat_id], a_start + a_duration)
                data_gen = self.data_gen_func(s.times[sat_id], a_start + a_duration)
                if (
                    s.energy[sat_id] + energy_gen + self.payload_power_gen * a_duration
                    < self.min_energy
                ):
                    continue
                if (
                    s.data[sat_id] + data_gen + self.payload_data_gen * a_duration
                    > self.max_data
                ):
                    continue

                # Consider action
                actions.append(
                    Action(
                        satellite_id=sat_id,
                        start=a_start,
                        duration=a_duration,
                        request=req,
                    )
                )
                actions_count += 1

                if actions_count >= N_max:
                    break

        return actions

    def get_discrete_index(self, ts, te):
        """
        Get the discrete index for a given time range

        Args:
            ts (float): Start time
            te (float): End time

        Returns:
            tuple[int, int]: Discrete index
        """

        i_s = int(ts / self.time_step)
        i_e = int(te / self.time_step)
        if i_e == ts:
            i_e += 1
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

        # No request
        if a.request is None:
            sp = deepcopy(s)
            sp.times[sat_id] = time
            return sp

        user_id = a.request.user_id
        request_id = a.request.id

        duration = min(
            a.duration,
            a.request.duration - s.request_times[request_id],
        )
        payload_on = user_id >= 0
        data = max(
            self.min_data,
            s.data[sat_id]
            + self.payload_data_gen * duration * payload_on
            + self.data_gen_func(s.times[sat_id], time),
        )
        energy = min(
            self.max_energy,
            s.energy[sat_id]
            + self.payload_power_gen * duration * payload_on
            + self.energy_gen_func(s.times[sat_id], time),
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
        sp.request_times[request_id] = s.request_times[request_id] + duration
        assert sp.request_times[request_id] <= a.request.duration
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
        return np.sum(cn0) * self.time_step

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

        payload_on = a.request.user_id >= 0
        if payload_on:
            return self.integrate_normalized_CN0(
                a.start, a.start + a.duration, a.request.user_id, a.satellite_id
            )
        else:
            return 0

    def initial_state(self) -> State:
        """
        Get initial state

        Returns:
            State: Initial state
        """

        return State(
            times=[0.0 for _ in range(self.N_satellites)],
            requests=[None for _ in range(self.N_satellites)],
            request_times={req.id: 0.0 for req in self.requests},
            data=[
                (self.min_data + self.max_data) / 2.0 for _ in range(self.N_satellites)
            ],
            energy=[
                (self.min_energy + self.max_energy) / 2.0
                for _ in range(self.N_satellites)
            ],
        )

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
        discounts = np.array([gamma**i for i in range(len(policy[:-1]))])
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
            sf.request_times[req.id] / req.duration * 100
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
