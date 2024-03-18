import numpy as np
import logging
from dataclasses import dataclass, field
from copy import deepcopy


@dataclass(frozen=True, repr=True)
class Request:
    id: int  # Request id
    user_id: int  # User id
    rv: np.array = field(repr=False)  # Position and velocity
    start: float  # Start time
    end: float  # End time
    duration: float  # Duration
    priority: int = 1  # Priority


@dataclass(frozen=True, repr=True)
class ServiceWindow:
    id: int  # ServiceWindow id
    satellite_id: int  # Satellite id
    request_id: int  # Request
    start: float  # Start time
    end: float  # End time


@dataclass(frozen=True, repr=False)
class State:
    time: list[float]  # Current time
    last_window: list[ServiceWindow]  # Last window
    request_time: dict[int, float]  # Request time
    data: list[float]  # Data onboard
    energy: list[float]  # Energy onboard

    def __repr__(self):
        s = "State("
        s += f"t={np.round(self.time, 2)}, "
        s += f"d={np.round(self.data, 2)}, "
        s += f"e={np.round(self.energy, 2)}, "
        s += f"request_time={self.request_time}"
        s += ")"
        return s

    def __hash__(self):
        return hash(
            (
                self.time,
                self.last_window,
                tuple(self.request_time.items()),
                self.data,
                self.energy,
            )
        )


@dataclass(frozen=True, repr=False)
class Action:
    satellite_id: int  # Satellite id
    start: float  # Start time
    duration: float  # Duration
    window: ServiceWindow  # ServiceWindow

    def __repr__(self):
        a = "Action("
        a += f"sat={self.satellite_id}, "
        a += f"request={self.window.request_id}, "
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
        # Hyperparameters
        N_max_actions: int,  # Maximum number of actions to consider
        min_action_duration: float,  # Minimum action duration
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
        self.request_dict = {r.id: r for r in requests}
        self.service_windows = sorted(service_windows, key=lambda w: w.start)
        self.transition_times = transition_times
        self.CN0_norm = CN0 / np.nanmax(CN0, axis=1)[:, None]
        self.N_satellites = len(set(w.satellite_id for w in service_windows))

        self.N_max_actions = N_max_actions
        self.min_action_duration = min_action_duration

        self.min_energy = min_energy
        self.max_energy = max_energy
        self.min_data = min_data
        self.max_data = max_data
        self.payload_power_gen = payload_power_gen
        self.payload_data_gen = payload_data_gen
        self.energy_gen_func = energy_gen_func
        self.data_gen_func = data_gen_func
        self.tf = max([r.end for r in requests])

    def available_actions(self, s: State) -> list[Action]:
        # Satellite to select an action
        sat_id = np.argmin(s.time)

        # List of available windows
        other_sat_request_ids = [
            w.request_id
            for other_sat_id, w in enumerate(s.last_window)
            if w is not None and other_sat_id != sat_id
        ]
        windows = [
            w
            for w in self.service_windows
            if
            # There is still time to start the window
            w.end >= s.time[sat_id]
            # The request is not completed
            and s.request_time[w.request_id] <= self.request_dict[w.request_id].duration
            # The window is for the selected satellite
            and w.satellite_id == sat_id
            # The request is not being served by another satellite
            and w.request_id not in other_sat_request_ids
        ]

        actions = []
        actions_count = 0
        for window in windows:
            # Time
            trans_time = (
                0
                if s.last_window[sat_id] is None
                else self.transition_times[
                    sat_id, s.last_window[sat_id].request_id, window.request_id
                ]
            )

            a_start = max(
                window.start,
                s.time[sat_id] + trans_time,
            )
            a_duration = min(
                self.min_action_duration,
                window.end - a_start,  # Window end
                self.request_dict[window.request_id].duration
                - s.request_time[window.request_id],  # Requested duration
            )
            if a_duration <= 0:
                continue

            # Resources
            energy_gen = self.energy_gen_func(s.time[sat_id], a_start + a_duration)
            data_gen = self.data_gen_func(s.time[sat_id], a_start + a_duration)
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
                    window=window,
                )
            )
            actions_count += 1

            if actions_count >= self.N_max_actions:
                break

        return actions

    def get_discrete_index(self, ts, te):
        # Get the discrete index for a given time range
        i_s = int(ts / self.time_step)
        i_e = int(te / self.time_step)
        if i_e == ts:
            i_e += 1
        return i_s, i_e

    def transition_function(self, s: State, a: Action) -> State:
        window = a.window
        time = a.start + a.duration
        sat_id = a.satellite_id

        payload_on = window.request_id >= 0
        data = max(
            self.min_data,
            s.data[sat_id]
            + self.payload_data_gen * a.duration * payload_on
            + self.data_gen_func(s.time[sat_id], time),
        )
        energy = min(
            self.max_energy,
            s.energy[sat_id]
            + self.payload_power_gen * a.duration * payload_on
            + self.energy_gen_func(s.time[sat_id], time),
        )

        eps = 1e-6
        assert s.time[sat_id] <= a.start + eps
        assert data >= self.min_data - eps
        assert data <= self.max_data + eps
        assert energy <= self.max_energy + eps
        assert energy >= self.min_energy - eps

        new_state = deepcopy(s)
        new_state.time[sat_id] = time
        new_state.last_window[sat_id] = window
        new_state.request_time[window.request_id] = min(
            s.request_time[window.request_id] + a.duration,
            self.request_dict[window.request_id].duration,
        )
        new_state.data[sat_id] = data
        new_state.energy[sat_id] = energy
        return new_state

    def old_reward_function(self, s: State, a: Action) -> float:
        return (
            a.window.reward
            * a.duration
            / self.request_dict[a.window.request_id].duration
        )

    def integrate_normalized_CN0(self, ts, te, request_id, satellite_id):
        # Integrate normalized CN0
        i_s, i_e = self.get_discrete_index(ts, te)
        cn0 = self.CN0_norm[satellite_id, request_id, i_s:i_e]
        cn0[np.isnan(cn0)] = 0
        return np.sum(cn0) * self.time_step

    def reward_function(self, s: State, a: Action) -> float:
        payload_on = a.window.request_id >= 0
        if payload_on:
            return self.integrate_normalized_CN0(
                a.start, a.start + a.duration, a.window.request_id, a.satellite_id
            )
        else:
            return 0

    def initial_state(self):
        return State(
            time=[0.0 for _ in range(self.N_satellites)],
            last_window=[None for _ in range(self.N_satellites)],
            request_time={r.id: 0.0 for r in self.request_dict.values()},
            data=[
                (self.min_data + self.max_data) / 2.0 for _ in range(self.N_satellites)
            ],
            energy=[
                (self.min_energy + self.max_energy) / 2.0
                for _ in range(self.N_satellites)
            ],
        )

    def total_reward(self, policy: list[tuple[State, Action]]) -> float:
        return sum(self.reward_function(s, a) for s, a in policy if a is not None)

    def percentage_completed(self, policy: list[tuple[State, Action]]) -> float:
        s0 = policy[0][0]
        sf = policy[-1][0]
        return {
            r.id: round(sf.request_time[r.id] / r.duration * 100, 2)
            for r in self.request_dict.values()
            if r.id >= 0
        }
