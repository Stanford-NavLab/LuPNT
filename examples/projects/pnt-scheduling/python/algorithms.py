import numpy as np
import cvxpy as cp
import logging
from dataclasses import dataclass, field
import utils


@dataclass(frozen=True, repr=True)
class Request:
    id: int  # Request id
    user_id: int  # User id
    rv: np.array = field(repr=False)  # Position and velocity
    start: float  # Start time
    end: float  # End time
    duration: float  # Duration


@dataclass(frozen=True, repr=True)
class ServiceWindow:
    id: int  # ServiceWindow id
    request_id: int  # Request
    start: float  # Start time
    end: float  # End time
    reward: float  # Reward
    power_gen: float  # Power generation
    data_gen: float  # Data generation


@dataclass(frozen=True, repr=False)
class State:
    time: float  # Current time
    last_window: ServiceWindow  # Last window
    request_time: dict[int, float]  # Request time
    data: float  # Data onboard
    energy: float  # Energy onboard

    def __repr__(self):
        s = "State("
        s += f"t={self.time:.2f}, "
        s += f"d={self.data:.2f}, "
        s += f"e={self.energy:.2f}, "
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
    start: float  # Start time
    duration: float  # Duration
    window: ServiceWindow  # ServiceWindow

    def __repr__(self):
        return f"Action(window={self.window.id}, start={self.start:.2f}, duration={self.duration:.2f})"


def logging_decorator(func):
    def wrapper(*args, **kwargs):
        logging.debug(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper


@dataclass()
class PntSchedulingProblem:

    def __init__(
        self,
        requests: list[Request],  # List of requests
        service_windows: list[ServiceWindow],  # List of service windows
        transition_times: np.ndarray,  # Transition times matrix
        N_max_actions: int,  # Maximum number of actions to consider
        min_action_duration: float,  # Minimum action duration
        min_energy: float,  # Minimum energy
        max_energy: float,  # Maximum energy
        min_data: float,  # Minimum data
        max_data: float,  # Maximum data
    ):
        self.request_dict = {r.id: r for r in requests}
        self.service_windows = sorted(service_windows, key=lambda w: w.start)
        self.transition_times = transition_times
        self.N_max_actions = N_max_actions
        self.min_action_duration = min_action_duration
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.min_data = min_data
        self.max_data = max_data

    def available_actions(self, s: State) -> list[ServiceWindow]:
        windows = [
            w
            for w in self.service_windows
            if w.end >= s.time  # There is still time to start the window
            and s.request_time[w.request_id]
            <= self.request_dict[w.request_id].duration  # The request is not completed
        ]

        actions = []
        actions_count = 0
        for window in windows:
            # Time
            trans_time = (
                0
                if s.last_window is None
                else self.transition_times[s.last_window.id, window.id]
            )
            action_start = max(window.start, s.time + trans_time)
            action_duration = min(
                self.min_action_duration,
                window.end - action_start,
                self.request_dict[window.request_id].duration
                - s.request_time[window.request_id],
            )
            if action_duration <= 0:
                continue

            # Resources

            # Consider action
            actions.append(
                Action(start=action_start, duration=action_duration, window=window)
            )
            actions_count += 1

            if actions_count >= self.N_max_actions:
                break

        return actions

    def transition_function(self, s: State, a: Action) -> State:
        window = a.window
        request_time = s.request_time.copy()
        request_time[window.request_id] += a.duration
        return State(
            time=a.start + a.duration,
            last_window=window,
            request_time=request_time,
            data=s.data + window.data_gen,
            energy=s.energy + window.power_gen,
        )

    def reward_function(self, s: State, a: Action) -> float:
        return (
            a.window.reward
            * a.duration
            / self.request_dict[a.window.request_id].duration
        )

    def initial_state(self):
        return State(
            time=0,
            last_window=None,
            request_time={r.id: 0 for r in self.request_dict.values()},
            data=0,
            energy=0,
        )

    def total_reward(self, policy: list[tuple[State, Action]]) -> float:
        return sum(self.reward_function(s, a) for s, a in policy if a is not None)


class SmdpForwardSearchSolver:
    """
    Forward search for solving Semi-Markov Decision Processes

    Args:
        problem (PntSchedulingProblem): Problem instance
    """

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

    def select_action(
        self, s: State, d: int, gamma: float
    ) -> tuple[ServiceWindow, float]:
        actions = self.problem.available_actions(s)
        if d == 0 or not actions:
            return None, 0
        a_star, v_star = None, -np.inf
        for a in actions:
            sp = self.problem.transition_function(s, a)
            _, vp = self.select_action(sp, d - 1, gamma)
            v = self.problem.reward_function(s, a) + gamma ** (a.start - s.time) * vp
            if v > v_star:
                a_star, v_star = a, v
        return a_star, v_star

    def solve(self, s: State, d: int, gamma: float) -> list[tuple[State, Action]]:
        logging.debug(f"Initial state: {s}")
        policy = []
        a, _ = self.select_action(s, d, gamma)
        while a is not None:
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
            a, _ = self.select_action(s, d, gamma)
        policy.append((s, None))
        return policy


class SmdpMctsSolver:
    """
    Monte Carlo Tree Search for solving Semi-Markov Decision Processes

    Args:
        problem (PntSchedulingProblem): Problem instance
    """

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

        self.N = dict[State, dict[Action, float]]()
        self.Q = dict[State, dict[Action, float]]()

    def get_N(self, s: State, a: Action) -> float:
        if s not in self.N or a not in self.N[s]:
            return 1
        return self.N[s][a]

    def get_Q(self, s: State, a: Action) -> float:
        if s not in self.Q or a not in self.Q[s]:
            return self.problem.reward_function(s, a)
        return self.Q[s][a]

    def rollout(self, s: State, d: int, gamma: float) -> float:
        actions = self.problem.available_actions(s)

        if d == 0 or not actions:
            return 0

        use_rewards = True
        if use_rewards:
            rewards = np.array([self.problem.reward_function(s, a) for a in actions])
            rewards += 1e-5  # Add small value to avoid division by zero
            rewards /= np.sum(rewards)
            idx = np.random.choice(len(actions), p=rewards)
        else:
            idx = np.random.choice(len(actions))

        a = actions[idx]
        r = self.problem.reward_function(s, a)
        sp = self.problem.transition_function(s, a)
        return r + gamma ** (a.start - s.time) * self.rollout(sp, d - 1, gamma)

    def bonus(self, c: float, N_s: int, N_sa: int) -> float:
        return c * np.sqrt(np.log(N_s) / N_sa) if N_sa > 0 else np.inf

    def simulate(self, s: State, d: int, gamma: float, c: float) -> float:
        actions = self.problem.available_actions(s)

        if d == 0 or not actions:
            return 0

        if s not in self.N:
            self.N[s] = dict[Action, float]()
            self.Q[s] = dict[Action, float]()
            for a in actions:
                self.N[s][a] = self.get_N(s, a)
                self.Q[s][a] = self.get_Q(s, a)
            return self.rollout(s, d, gamma)

        N_s = sum(self.N[s][a] for a in actions)
        a = max(actions, key=lambda a: self.Q[s][a] + self.bonus(c, N_s, self.N[s][a]))
        sp = self.problem.transition_function(s, a)
        r = self.problem.reward_function(s, a)
        q = r + gamma ** (a.start - s.time) * self.simulate(sp, d - 1, gamma, c)
        self.N[s][a] += 1
        self.Q[s][a] += (q - self.Q[s][a]) / self.N[s][a]
        return q

    def solve(
        self, s: State, d: int, gamma: float, n: int, c: float
    ) -> list[tuple[State, Action]]:
        """
        Solve the problem

        Args:
            s (State): Initial state
            d (int): Depth
            gamma (float): Discount factor
            n (int): Number of simulations
            c (float): Exploration constant
        """
        np.random.seed(0)

        for _ in range(n):
            self.simulate(s, d, gamma, c)

        policy = []
        actions = self.problem.available_actions(s)
        while actions:
            a = max(actions, key=lambda a: self.get_Q(s, a))
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
            actions = self.problem.available_actions(s)

        policy.append((s, None))
        return policy


class DiscreteTimeIpSolver:
    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

    def solve(self, state: State, time_step: float) -> list[tuple[State, Action]]:

        N_requests = len(self.problem.request_dict)
        N_time_steps = int(
            max([r.end for r in self.problem.request_dict.values()]) / time_step
        )
        durations = np.array(
            [r.duration / time_step for r in self.problem.request_dict.values()]
        )
        transition_times = np.ceil(self.problem.transition_times / time_step).astype(
            int
        )
        rewards = np.zeros((N_requests, N_time_steps))
        service_windows = np.zeros((N_requests, N_time_steps), dtype=bool)
        for i, r in enumerate(self.problem.request_dict.values()):
            windows = [w for w in self.problem.service_windows if w.request_id == r.id]
            starts = np.ceil([w.start / time_step for w in windows]).astype(int)
            ends = np.floor([w.end / time_step for w in windows]).astype(int)
            for w, s, e in zip(windows, starts, ends):
                service_windows[i, s:e] = True
                rewards[i, s:e] = w.reward / (r.duration / time_step)

        # Variables
        x = cp.Variable((N_requests, N_time_steps), boolean=True)

        # Objective
        objective = cp.Maximize(cp.sum(cp.multiply(x, rewards)))

        # Constraints
        constraints = []
        constraints.append(cp.sum(x, axis=0) <= 1)  # One action at a time
        constraints.append(cp.sum(x, axis=1) <= durations)  # Duration
        constraints.append(x <= service_windows)  # Service window
        for i in range(N_requests):
            for j in range(N_requests):
                if i == j:
                    continue
                for t in range(N_time_steps - 1):
                    tt = transition_times[i, j]
                    constraints.append(
                        x[j, t : t + tt + 1] <= 1 - x[i, t]
                    )  # Transition times

        # Optimize
        problem = cp.Problem(objective, constraints)
        problem.solve()
        print(problem.status)

        # Policy
        policy = self.construct_policy(state, x, time_step)
        return policy

    def construct_policy(
        self, state: State, x: np.array, time_step: float
    ) -> list[tuple[State, Action]]:

        solution = np.round(x.value).astype(int)  # Convert to integer
        service_windows_dict: dict[int, list[ServiceWindow]] = {
            r.id: [w for w in self.problem.service_windows if w.request_id == r.id]
            for r in self.problem.request_dict.values()
        }

        def get_window(s: int, e: int, request_id: int) -> int:
            windows = [
                w
                for w in service_windows_dict[request_id]
                if s >= int(np.ceil(w.start / time_step))
                and e <= int(np.floor(w.end / time_step))
            ]
            return windows[0]

        actions: list[Action] = []
        for i, r in enumerate(self.problem.request_dict.values()):
            start, end = utils.get_start_end_indexes(solution[i])
            for s, e in zip(start, end):
                window = get_window(s, e, r.id)
                a = Action(
                    start=s * time_step, duration=(e - s) * time_step, window=window
                )
                actions.append(a)
        actions.sort(key=lambda a: a.start)

        policy = []
        s = state
        for a in actions:
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
        policy.append((s, None))
        return policy
