import numpy as np
from tqdm.notebook import tqdm
from problem import PntSchedulingProblem, State, Action, ServiceWindow
from .Solver import Solver
from . import RuleBasedSolver

import networkx as nx


class SmdpForwardSearchSolver(Solver):
    """
    Forward search for solving Semi-Markov Decision Processes

    Args:
        problem (PntSchedulingProblem): Problem instance
        keep_tree (bool): Keep the search tree
    """

    def __init__(self, problem: PntSchedulingProblem, keep_tree: bool = True):
        self.problem = problem

        self.gamma = None
        self.N_max = None
        self.dur_min = None
        self.use_rule_based = False
        self.tree = nx.DiGraph() if keep_tree else None

    def select_action(self, s: State, d: int) -> tuple[ServiceWindow, float]:
        actions = self.problem.available_actions(s, self.N_max, self.dur_min)
        if d == 0 or not actions:
            if self.use_rule_based:
                policy = RuleBasedSolver(self.problem).solve(
                    s, N_max=self.N_max, dur_min=self.dur_min
                )
                return None, self.problem.total_reward(policy, gamma=self.gamma)
            return None, 0

        a_star, v_star = None, -np.inf

        if self.tree is not None:
            self.tree.add_node(s, value=v_star)

        for a in actions:
            sp = self.problem.transition_function(s, a)
            sat_id = a.sat_id
            _, vp = self.select_action(sp, d - 1)
            r = self.problem.reward_function(s, a)
            v = r + self.gamma ** (a.ts - s.t[sat_id]) * vp
            if v > v_star:
                a_star, v_star = a, v

            if self.tree is not None:
                self.tree.add_node(sp, value=v)
                self.tree.add_node(a, reward=r, value=vp)
                self.tree.add_edge(s, a)
                self.tree.add_edge(a, sp)

        self.tree.nodes[s]["value"] = v_star

        return a_star, v_star

    def solve(
        self,
        s: State,
        d: int,
        gamma: float,
        N_max: int,
        dur_min: float,
        use_rule_based: bool = False,
        progress: bool = False,
    ) -> list[tuple[State, Action]]:
        """
        Solve the problem using forward search

        Args:
            s (State): Initial state
            d (int): Depth of the search
            gamma (float): Discount factor
            N_max (int): Maximum number of actions
            dur_min (float): Minimum action duration
            progress (bool): Show progress bar

        Returns:
            list[tuple[State, Action]]: Policy
        """

        self.gamma = gamma
        self.N_max = N_max
        self.dur_min = dur_min
        self.use_rule_based = use_rule_based

        policy = []
        a, _ = self.select_action(s, d)
        policy.append((s, a))

        # Progress bar
        if progress:
            tf = self.problem.t_final
            t = min(s.t)
            bar = tqdm(
                total=int(tf - t), desc="Solving Forward Search (progress in hours)"
            )

        while a is not None:
            s = self.problem.transition_function(s, a)
            a, _ = self.select_action(s, d)
            policy.append((s, a))

            # Update progress bar
            if progress:
                bar.update(int(min(s.t) - t))
                t = min(s.t)

        # Close progress bar
        if progress:
            bar.update(int(tf - bar.n))

        policy = self.problem.clean_policy(policy)
        return policy