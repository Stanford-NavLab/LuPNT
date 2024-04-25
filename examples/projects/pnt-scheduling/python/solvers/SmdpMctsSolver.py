import numpy as np
from tqdm.notebook import tqdm
from problem import PntSchedulingProblem, State, Action
from .Solver import Solver


class SmdpMctsSolver(Solver):
    """
    Monte Carlo Tree Search for solving Semi-Markov Decision Processes

    Args:
        problem (PntSchedulingProblem): Problem instance
    """

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

        self.N = dict[State, dict[Action, int]]()
        self.Q = dict[State, dict[Action, float]]()

        self.c = None
        self.gamma = None
        self.N_max = None
        self.d_min = None

    def get_N(self, s: State, a: Action) -> int:
        if s not in self.N or a not in self.N[s]:
            return 0
        return self.N[s][a]

    def get_Q(self, s: State, a: Action) -> float:
        if s not in self.Q or a not in self.Q[s]:
            return self.problem.reward_function(s, a)
        return self.Q[s][a]

    def rollout(self, s: State, d: int) -> float:
        actions = self.problem.available_actions(s, self.N_max, self.d_min)

        if d == 0 or not actions:
            return 0

        use_rewards = True
        if use_rewards:
            rewards = np.array([self.problem.reward_function(s, a) for a in actions])
            rewards = rewards + 1e-5  # Add small value to avoid division by zero
            rewards /= np.sum(rewards)
            idx = np.random.choice(len(actions), p=rewards)
        else:
            idx = np.random.choice(len(actions))

        a = actions[idx]
        sat_id = a.sat_id
        r = self.problem.reward_function(s, a)
        sp = self.problem.transition_function(s, a)
        return r + self.gamma ** (a.ts - s.t[sat_id]) * self.rollout(sp, d - 1)

    def bonus(self, N_s: int, N_sa: int) -> float:
        return self.c * np.sqrt(np.log(N_s) / N_sa) if N_sa > 0 else np.inf

    def simulate(self, s: State, d: int) -> float:
        actions = self.problem.available_actions(s, self.N_max, self.d_min)

        if d == 0 or not actions:
            return 0

        if s not in self.N:
            self.N[s] = dict[Action, float]()
            self.Q[s] = dict[Action, float]()
            for a in actions:
                self.N[s][a] = self.get_N(s, a)
                self.Q[s][a] = self.get_Q(s, a)
            return self.rollout(s, d)

        N_s = sum(self.N[s][a] for a in actions)
        a = max(actions, key=lambda a: self.Q[s][a] + self.bonus(N_s, self.N[s][a]))
        sat_id = a.sat_id

        sp = self.problem.transition_function(s, a)
        r = self.problem.reward_function(s, a)
        q = r + self.gamma ** (a.ts - s.t[sat_id]) * self.simulate(sp, d - 1)

        self.N[s][a] += 1
        self.Q[s][a] += (q - self.Q[s][a]) / self.N[s][a]
        return q

    def solve(
        self,
        s: State,
        d: int,
        gamma: float,
        n: int,
        c: float,
        N_max: int,
        d_min: float,
        progress: bool = False,
    ) -> list[tuple[State, Action]]:
        """
        Solve the problem

        Args:
            s (State): Initial state
            d (int): Depth
            gamma (float): Discount factor
            n (int): Number of simulations
            c (float): Exploration constant
            N_max (int): Maximum number of actions
            d_min (float): Minimum action duration
            progress (bool): Show progress bar

        Returns:
            list[tuple[State, Action]]: Policy
        """

        self.c = c
        self.gamma = gamma
        self.N_max = N_max
        self.d_min = d_min

        for _ in range(n):
            self.simulate(s, d)

        # Progress bar
        if progress:
            tf = self.problem.t_final
            t = min(s.t)
            bar = tqdm(total=int(tf - t), desc="Solving MCTS (progress in hours)")

        policy = []
        actions = self.problem.available_actions(s, self.N_max, self.d_min)
        while actions:
            a = max(actions, key=lambda a: self.get_Q(s, a))
            policy.append((s, a))
            s = self.problem.transition_function(s, a)

            for _ in range(n):
                self.simulate(s, d)
            actions = self.problem.available_actions(s, self.N_max, self.d_min)

            if progress:
                # Update progress bar
                bar.update(int(min(s.t) - t))
                t = min(s.t)

        if progress:
            # Close progress bar
            bar.update(int(tf - bar.n))

        policy.append((s, None))
        policy = self.problem.clean_policy(policy)
        return policy
