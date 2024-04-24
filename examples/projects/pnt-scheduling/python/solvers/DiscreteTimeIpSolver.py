import numpy as np
import cvxpy as cp
import util
import logging
import gurobipy
from problem import PntSchedulingProblem, State, Action, ServiceWindow
from .Solver import Solver


class DiscreteTimeIpSolver(Solver):
    def __init__(self, problem: PntSchedulingProblem):
        self.problem: PntSchedulingProblem = problem
        self.solution: np.array = None

    def solve(self, s: State, time_step_factor: float) -> list[tuple[State, Action]]:

        N_sat = self.problem.N_sat
        N_req = len(self.problem.request_dict)
        Dt = self.problem.t_final / self.problem.CN0_norm.shape[2]
        dt = time_step_factor * Dt
        N_t = int(self.problem.t_final / dt)

        durations = np.ceil(
            [r.duration / dt for r in self.problem.request_dict.values()]
        ).astype(int)
        transition_times = np.ceil(self.problem.transition_times / dt).astype(int)
        rewards = np.zeros((N_sat, N_req, N_t))
        service_windows = np.zeros((N_sat, N_req, N_t))
        data_gen_requests = np.zeros((N_sat, N_req, N_t))
        energy_gen_requests = np.zeros((N_sat, N_req, N_t))
        for i_sat in range(N_sat):
            for j_req, req in enumerate(self.problem.request_dict.values()):
                windows = [
                    win
                    for win in self.problem.service_windows
                    if win.request_id == req.id and win.satellite_id == i_sat
                ]
                starts = np.ceil([w.start / dt for w in windows]).astype(int)
                ends = np.floor([w.end / dt for w in windows]).astype(int)
                for w, ts, te in zip(windows, starts, ends):
                    if ts >= te:
                        logging.debug(
                            f"Service window {w.id} has duration 0. Try decreasing time step."
                        )
                        continue
                        # if te < N_time_steps:
                        #     te += 1
                        # else:
                        #     ts -= 1
                    service_windows[i_sat, j_req, ts:te] = 1
                    payload_on = w.request_id >= 0
                    f = time_step_factor
                    if payload_on:
                        rewards[i_sat, j_req, ts:te] = (
                            self.problem.CN0_norm[
                                i_sat, w.request_id, (f * ts) : (f * te) : f
                            ]
                            * dt
                        )
                        assert not np.isnan(rewards).any()
                        data_gen_requests[i_sat, j_req, ts:te] = (
                            self.problem.payload_data_gen * dt
                        )
                        energy_gen_requests[i_sat, j_req, ts:te] = (
                            self.problem.payload_energy_gen * dt
                        )

        data_gen = self.problem.data_gen_func(
            np.arange(N_t) * dt,
            (np.arange(N_t) + 1) * dt,
        )
        energy_gen = self.problem.energy_gen_func(
            np.arange(N_t) * dt,
            (np.arange(N_t) + 1) * dt,
        )
        # Variables
        x = [cp.Variable((N_req, N_t), boolean=True) for _ in range(N_sat)]

        # Objective
        lambda_sep = 1e-5 / (N_req * N_t)
        inner_sum = [
            cp.sum(
                cp.multiply(x[i], rewards[i])
                - lambda_sep * cp.sum(cp.abs(cp.diff(x[i], axis=1)))
            )
            for i in range(N_sat)
        ]
        objective = cp.Maximize(cp.sum(inner_sum))

        # Constraints
        constraints = []
        for i_sat in range(N_sat):
            tmp_data = data_gen + cp.sum(
                cp.multiply(data_gen_requests[i_sat], x[i_sat]), axis=0
            )
            tmp_energy = energy_gen + cp.sum(
                cp.multiply(energy_gen_requests[i_sat], x[i_sat]), axis=0
            )
            for t in range(N_t + 1):
                constraints.append(
                    s.data[i_sat] + cp.sum(tmp_data[:t]) <= self.problem.max_data
                )
                constraints.append(
                    s.energy[i_sat] + cp.sum(tmp_energy[:t]) >= self.problem.min_energy
                )

        # One action at a time for each satellite
        constraints.extend([cp.sum(x[i_sat], axis=0) <= 1 for i_sat in range(N_sat)])

        # One satellite at a time for each request
        constraints.extend(
            [
                cp.sum([x[i_sat][j_req, t] for i_sat in range(N_sat)]) <= 1
                for t in range(N_t)
                for j_req, req in enumerate(self.problem.request_dict.values())
                if req.id >= 0
            ]
        )

        # Duration
        inner_sum = [cp.sum(x[i_sat], axis=1) for i_sat in range(N_sat)]
        constraints.append(cp.sum(inner_sum) <= durations)

        # Service window
        constraints.extend(
            [x[i_sat] <= service_windows[i_sat] for i_sat in range(N_sat)]
        )
        for i_sat in range(N_sat):
            for i in range(N_req):
                for j in range(N_req):
                    if i == j:
                        continue
                    for t in range(N_t - 1):
                        tt = transition_times[i_sat, i, j]
                        # Transition times
                        constraints.append(
                            x[i_sat][j, t : t + tt + 1] <= 1 - x[i_sat][i, t]
                        )

        # Optimize
        env = gurobipy.Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        problem = cp.Problem(objective, constraints)
        # verbose if logging is set to debug
        problem.solve(
            solver=cp.GUROBI, env=env, verbose=logging.DEBUG <= logging.root.level
        )

        # Convert to integer
        self.solution = np.zeros((N_sat, N_req, N_t))
        for i_sat in range(N_sat):
            self.solution[i_sat] = np.round(x[i_sat].value).astype(int)

        # Policy
        # policy = self.construct_policy(s, dt)
        # return policy

        # def construct_policy(
        #     self, s: State, time_step: float
        # ) -> list[tuple[State, Action]]:

        # Dictionary mapping request and satellite ids to service windows
        service_windows_dict: dict[int, dict[int, list[ServiceWindow]]] = {
            req.id: {
                i_sat: [
                    w
                    for w in self.problem.service_windows
                    if w.request_id == req.id and w.satellite_id == i_sat
                ]
                for i_sat in range(self.problem.N_sat)
            }
            for req in self.problem.request_dict.values()
        }

        # Get window id given a time interval and a request id
        def get_window(ts: int, te: int, request_id: int, satellite_id: int) -> int:
            windows = [
                win
                for win in service_windows_dict[request_id][satellite_id]
                if ts >= int(np.ceil(win.start / dt)) - 1
                and te <= int(np.floor(win.end / dt)) + 1
            ]
            return windows[0]

        # Create actions by iterating over the requests
        actions: list[Action] = []
        for j_req, req in enumerate(self.problem.request_dict.values()):
            for i_sat in range(self.problem.N_sat):
                start, end = util.get_start_end_indexes(self.solution[i_sat, j_req])
                for ts, te in zip(start, end):
                    window = get_window(ts, te, req.id, i_sat)
                    a = Action(
                        satellite_id=i_sat,
                        start=ts * dt,
                        duration=(te - ts) * dt,
                        window=window,
                    )
                    actions.append(a)
        actions.sort(key=lambda a: a.start)

        tmp_d = [
            data_gen
            + np.sum(
                np.multiply(data_gen_requests[i_sat], self.solution[i_sat]), axis=0
            )
            for i_sat in range(N_sat)
        ]
        tmp_e = [
            energy_gen
            + np.sum(
                np.multiply(energy_gen_requests[i_sat], self.solution[i_sat]), axis=0
            )
            for i_sat in range(N_sat)
        ]

        # Create policy
        policy = []
        for a in actions:
            # Fix the duration of the action
            policy.append((s, a))
            s = self.problem.transition_function(s, a)

        policy.append((s, None))
        return policy
