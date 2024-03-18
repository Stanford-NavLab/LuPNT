import numpy as np
import cvxpy as cp
import logging
from dataclasses import dataclass, field
import utils
from tqdm import tqdm
from problem import PntSchedulingProblem, State, Action, ServiceWindow


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
