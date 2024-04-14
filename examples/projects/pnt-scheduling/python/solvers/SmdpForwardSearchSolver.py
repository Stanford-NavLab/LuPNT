import numpy as np
from tqdm.notebook import tqdm
from problem import PntSchedulingProblem, State, Action, ServiceWindow
from .Solver import Solver


class SmdpForwardSearchSolver(Solver):
    """
    Forward search for solving Semi-Markov Decision Processes

    Args:
        problem (PntSchedulingProblem): Problem instance
    """

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

        self.gamma = None
        self.N_max = None
        self.d_min = None

    def select_action(self, s: State, d: int) -> tuple[ServiceWindow, float]:
        actions = self.problem.available_actions(s, self.N_max, self.d_min)
        if d == 0 or not actions:
            return None, 0
        a_star, v_star = None, -np.inf
        for a in actions:
            sp = self.problem.transition_function(s, a)
            sat_id = a.satellite_id
            _, vp = self.select_action(sp, d - 1)
            v = (
                self.problem.reward_function(s, a)
                + self.gamma ** (a.duration + 1e-3) * vp
            )
            if v > v_star:
                a_star, v_star = a, v
        return a_star, v_star

    def solve(
        self,
        s: State,
        d: int,
        gamma: float,
        N_max: int,
        d_min: float,
        progress: bool = False,
    ) -> list[tuple[State, Action]]:
        """
        Solve the problem using forward search

        Args:
            s (State): Initial state
            d (int): Depth of the search
            gamma (float): Discount factor
            N_max (int): Maximum number of actions
            d_min (float): Minimum action duration
            progress (bool): Show progress bar

        Returns:
            list[tuple[State, Action]]: Policy
        """

        self.gamma = gamma
        self.N_max = N_max
        self.d_min = d_min

        policy = []
        a, _ = self.select_action(s, d)
        policy.append((s, a))

        # Progress bar
        if progress:
            tf = self.problem.tf
            t = min(s.times)
            bar = tqdm(
                total=int(tf - t), desc="Solving Forward Search (progress in hours)"
            )

        while a is not None:
            s = self.problem.transition_function(s, a)
            a, _ = self.select_action(s, d)
            policy.append((s, a))

            # Update progress bar
            if progress:
                bar.update(int(min(s.times) - t))
                t = min(s.times)

        # Close progress bar
        if progress:
            bar.update(int(tf - bar.n))

        policy = self.problem.clean_policy(policy)
        return policy
