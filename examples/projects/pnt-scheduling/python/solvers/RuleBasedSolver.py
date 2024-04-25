import numpy as np
from problem import PntSchedulingProblem, State, Action, ServiceWindow
from .Solver import Solver


class RuleBasedSolver(Solver):

    def __init__(self, problem: PntSchedulingProblem):
        self.problem = problem

    def solve(self, s: State, N_max, dur_min) -> list[tuple[State, Action]]:
        """
        Solve the problem using a rule-based approach

        Args:
            s (State): Initial state

        Returns:
            list[tuple[State, Action]]: Policy
        """

        def sorting_key(a: Action) -> float:
            return a.ts * self.problem.t_final + a.req.te if a.req else np.inf

        actions = sorted(
            self.problem.available_actions(s, N_max=N_max, dur_min=dur_min),
            key=sorting_key,
        )
        policy = []
        while actions:
            a = actions[0]
            policy.append((s, a))
            s = self.problem.transition_function(s, a)
            actions = sorted(
                self.problem.available_actions(s, N_max=N_max, dur_min=dur_min),
                key=sorting_key,
            )
        policy.append((s, None))
        policy = self.problem.clean_policy(policy)
        return policy
