import numpy as np
import cvxpy as cp
import logging
from dataclasses import dataclass, field
import utils
from tqdm import tqdm


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
        energy_gen_func: callable,  # Energy generation function
        data_gen_func: callable,  # Data generation function
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
        self.energy_gen_func = energy_gen_func
        self.data_gen_func = data_gen_func
        self.tf = max([r.end for r in requests])

    def available_actions(self, s: State) -> list[Action]:
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
            a_start = max(window.start, s.time + trans_time)
            a_duration = min(
                self.min_action_duration,
                window.end - a_start,
                self.request_dict[window.request_id].duration
                - s.request_time[window.request_id],
            )
            if a_duration <= 0:
                continue

            # Resources
            energy_gen = self.energy_gen_func(s.time, a_start + a_duration)
            data_gen = self.data_gen_func(s.time, a_start + a_duration)
            if s.energy + energy_gen + window.power_gen < self.min_energy:
                continue
            if s.data + data_gen + window.data_gen > self.max_data:
                continue

            # Consider action
            actions.append(Action(start=a_start, duration=a_duration, window=window))
            actions_count += 1

            if actions_count >= self.N_max_actions:
                break

        return actions

    def transition_function(self, s: State, a: Action) -> State:
        window = a.window
        time = a.start + a.duration
        request_time = s.request_time.copy()
        request_time[window.request_id] = min(
            request_time[window.request_id] + a.duration,
            self.request_dict[window.request_id].duration,
        )
        data = max(
            self.min_data,
            s.data + window.data_gen * a.duration + self.data_gen_func(s.time, time),
        )
        energy = min(
            self.max_energy,
            s.energy
            + window.power_gen * a.duration
            + self.energy_gen_func(s.time, time),
        )

        eps = 1e-6
        assert s.time <= a.start + eps
        assert data >= self.min_data - eps
        assert data <= self.max_data + eps
        assert energy <= self.max_energy + eps
        assert energy >= self.min_energy - eps

        return State(
            time=time,
            last_window=window,
            request_time=request_time,
            data=data,
            energy=energy,
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
            data=(self.min_data + self.max_data) / 2,
            energy=(self.min_energy + self.max_energy) / 2,
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

        # Progress bar
        tf = self.problem.tf
        t = s.time
        bar = tqdm(total=int(tf - t), desc="Solving Forward Search (progress in hours)")

        while a is not None:
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
            a, _ = self.select_action(s, d, gamma)

            # Update progress bar
            bar.update(int(s.time - t))
            t = s.time

        # Close progress bar
        bar.update(int(tf - bar.n))

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

        self.N = dict[State, dict[Action, int]]()
        self.Q = dict[State, dict[Action, float]]()

    def get_N(self, s: State, a: Action) -> int:
        if s not in self.N or a not in self.N[s]:
            return 0
        return self.N[s][a]

    def get_Q(self, s: State, a: Action) -> float:
        if s not in self.Q or a not in self.Q[s]:
            return self.problem.reward_function(s, a)
        return self.Q[s][a]

    def rollout(self, s: State, d: int, gamma: float) -> float:
        actions = self.problem.available_actions(s)

        if d == 0 or not actions:
            return 0

        use_rewards = False
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

        # Progress bar
        tf = self.problem.tf
        t = s.time
        bar = tqdm(total=int(tf - t), desc="Solving Forward Search (progress in hours)")

        policy = []
        actions = self.problem.available_actions(s)
        while actions:
            a = max(actions, key=lambda a: self.get_Q(s, a))
            policy.append((s, a))
            s = self.problem.transition_function(s, a)

            for _ in range(n):
                self.simulate(s, d, gamma, c)
            actions = self.problem.available_actions(s)

            # Update progress bar
            bar.update(int(s.time - t))
            t = s.time

        # Close progress bar
        bar.update(int(tf - bar.n))

        policy.append((s, None))
        return policy


class DiscreteTimeIpSolver:
    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem
        self.solution: np.array = None

    def solve(self, s: State, time_step: float) -> list[tuple[State, Action]]:

        N_requests = len(self.problem.request_dict)
        N_time_steps = int(self.problem.tf / time_step)
        durations = np.ceil(
            [r.duration / time_step for r in self.problem.request_dict.values()]
        ).astype(int)
        transition_times = np.ceil(self.problem.transition_times / time_step).astype(
            int
        )
        rewards = np.zeros((N_requests, N_time_steps))
        service_windows = np.zeros((N_requests, N_time_steps))
        data_gen_requests = np.zeros((N_requests, N_time_steps))
        energy_gen_requests = np.zeros((N_requests, N_time_steps))
        for i, r in enumerate(self.problem.request_dict.values()):
            windows = [w for w in self.problem.service_windows if w.request_id == r.id]
            starts = np.ceil([w.start / time_step for w in windows]).astype(int)
            ends = np.floor([w.end / time_step for w in windows]).astype(int)
            for w, ts, te in zip(windows, starts, ends):
                if ts == te:
                    print(
                        f"Service window {w.id} has duration 0. Try decreasing time step."
                    )
                    # if te < N_time_steps:
                    #     te += 1
                    # else:
                    #     ts -= 1
                service_windows[i, ts:te] = 1
                rewards[i, ts:te] = w.reward / durations[i]
                data_gen_requests[i, ts:te] = w.data_gen * time_step
                energy_gen_requests[i, ts:te] = w.power_gen * time_step

        data_gen = self.problem.data_gen_func(
            np.arange(N_time_steps) * time_step,
            (np.arange(N_time_steps) + 1) * time_step,
        )
        energy_gen = self.problem.energy_gen_func(
            np.arange(N_time_steps) * time_step,
            (np.arange(N_time_steps) + 1) * time_step,
        )
        # Variables
        x = cp.Variable((N_requests, N_time_steps), boolean=True)

        # Objective
        lambda_sep = 1e-3 / (N_requests * N_time_steps)
        objective = cp.Maximize(
            cp.sum(
                cp.multiply(x, rewards)
                - lambda_sep * cp.sum(cp.abs(cp.diff(x, axis=1)))
            )
        )

        # Constraints
        constraints = []
        for t in range(N_time_steps):
            constraints.append(
                s.data
                + np.sum(data_gen[: t + 1])
                + cp.sum(cp.multiply(data_gen_requests, x)[:, : t + 1])
                <= self.problem.max_data
            )
            constraints.append(
                s.energy
                + np.sum(energy_gen[: t + 1])
                + cp.sum(cp.multiply(energy_gen_requests, x)[:, : t + 1])
                >= self.problem.min_energy
            )

        constraints.append(cp.sum(x, axis=0) <= 1)  # One action at a time
        constraints.append(cp.sum(x, axis=1) <= durations)  # Duration
        constraints.append(x <= service_windows)  # Service window
        for i in range(N_requests):
            for j in range(N_requests):
                if i == j:
                    continue
                for t in range(N_time_steps - 1):
                    tt = transition_times[i, j]
                    # Transition times
                    constraints.append(x[j, t : t + tt + 1] <= 1 - x[i, t])

        # Optimize
        problem = cp.Problem(objective, constraints)
        problem.solve(verbose=True, solver=cp.GUROBI)
        self.solution = np.round(x.value).astype(int)  # Convert to integer
        print(problem.status)

        # Policy
        policy = self.construct_policy(s, time_step)
        return policy

    def construct_policy(
        self, s: State, time_step: float
    ) -> list[tuple[State, Action]]:

        # Dictionary mapping request id to service windows
        service_windows_dict: dict[int, list[ServiceWindow]] = {
            r.id: [w for w in self.problem.service_windows if w.request_id == r.id]
            for r in self.problem.request_dict.values()
        }

        # Get window id given a time interval and a request id
        def get_window(ts: int, te: int, request_id: int) -> int:
            windows = [
                w
                for w in service_windows_dict[request_id]
                if ts >= int(np.ceil(w.start / time_step)) - 1
                and te <= int(np.floor(w.end / time_step)) + 1
            ]
            return windows[0]

        # Create actions by iterating over the requests
        actions: list[Action] = []
        for i, r in enumerate(self.problem.request_dict.values()):
            start, end = utils.get_start_end_indexes(self.solution[i])
            for ts, te in zip(start, end):
                window = get_window(ts, te, r.id)
                a = Action(
                    start=ts * time_step, duration=(te - ts) * time_step, window=window
                )
                actions.append(a)
        actions.sort(key=lambda a: a.start)

        # Create policy
        policy = []
        for a in actions:
            # Fix the duration of the action
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
        policy.append((s, None))
        return policy


class RuleBasedSolver:

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

    def solve(self, s: State) -> list[tuple[State, Action]]:
        def sorting_key(a: Action) -> float:
            return (
                a.start * self.problem.tf + a.window.request_id
                if a.window.request_id >= 0
                else np.inf
            )

        actions = sorted(self.problem.available_actions(s), key=sorting_key)
        policy = []
        while actions:
            a = actions[0]
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
            actions = sorted(self.problem.available_actions(s), key=sorting_key)
        policy.append((s, None))
        return policy
