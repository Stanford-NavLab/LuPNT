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
        s += f"requests={list(self.request_times.values())}"
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


@dataclass
class PntSchedulingProblem:

    # Scenario
    time_step: float  # Time step
    requests: list[Request]  # List of requests
    service_windows: list[ServiceWindow]  # List of service windows
    transition_times: np.ndarray  # Transition times matrix
    CN0: np.ndarray  # Carrier-to-noise density ratio [dB-Hz]

    # Resources
    min_energy: float  # Minimum energy
    max_energy: float  # Maximum energy
    min_data: float  # Minimum data
    max_data: float  # Maximum data
    payload_power_gen: float  # Payload power generation
    payload_data_gen: float  # Payload data generation
    energy_gen_func: callable  # Energy generation function
    data_gen_func: callable  # Data generation function

    # Online
    current_policy: list[tuple[State, Action]] = None  # Existing schedule
    current_time: float = 0  # Current time

    def __post_init__(self):
        self.service_windows_dict = {
            req.id: [win for win in self.service_windows if win.user_id == req.user_id]
            for req in self.requests
        }
        self.CN0_norm = self.CN0 / np.nanmax(self.CN0, axis=2)[:, :, None]
        self.N_satellites = len(set(w.satellite_id for w in self.service_windows))
        self.tf = max([r.end for r in self.requests])

    def set_current_time(self, current_time):
        self.current_time = current_time

    def set_current_policy(self, policy: list[tuple[State, Action]]):
        self.current_policy = deepcopy(policy)

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
            if req.arrival <= self.current_time
        }

        # List of available requests
        available_requests = [
            req
            for req in self.requests
            if req.arrival <= self.current_time
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
                duration=min(d_min, self.tf - s.times[sat_id]),
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
                    req.duration - s.request_times[req.id],  # Requested duration
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

        if a.request is not None:
            usr_id = a.request.user_id
            req_id = a.request.id
            payload_on = usr_id >= 0

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

        if self.current_policy is None:
            return State(
                times=[0.0 for _ in range(self.N_satellites)],
                requests=[None for _ in range(self.N_satellites)],
                request_times={req.id: 0.0 for req in self.requests},
                data=[
                    (self.min_data + self.max_data) / 2.0
                    for _ in range(self.N_satellites)
                ],
                energy=[
                    (self.min_energy + self.max_energy) / 2.0
                    for _ in range(self.N_satellites)
                ],
            )

        for s, a in self.current_policy:
            print(s.times, a.start, self.current_time)
            for sat_id in range(self.N_satellites):
                if (
                    s.times[sat_id] > self.current_time
                    or a is None
                    or a.start > self.current_time
                ):
                    return s
                pass

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

        actions_by_sat = [list[Action]() for _ in range(self.N_satellites)]
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
