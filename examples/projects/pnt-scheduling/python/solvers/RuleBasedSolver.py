import numpy as np
from problem import PntSchedulingProblem, State, Action, ServiceWindow
from .Solver import Solver


class RuleBasedSolver(Solver):

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

    def solve(self, s: State) -> list[tuple[State, Action]]:
        """
        Solve the problem using a rule-based approach

        Args:
            s (State): Initial state

        Returns:
            list[tuple[State, Action]]: Policy
        """

        def sorting_key(a: Action) -> float:
            return (
                a.start * self.problem.t_final + a.window.request_id
                if a.window.request_id >= 0
                else np.inf
            )

        N_max = 2
        d_min = 1

        actions = sorted(
            self.problem.available_actions(s, N_max=N_max, d_min=d_min), key=sorting_key
        )
        policy = []
        while actions:
            a = actions[0]
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
            actions = sorted(
                self.problem.available_actions(s, N_max=N_max, d_min=d_min),
                key=sorting_key,
            )
        policy.append((s, None))
        return policy
