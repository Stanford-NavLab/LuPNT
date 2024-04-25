import numpy as np
import cvxpy as cp
import util
import logging
import gurobipy
from problem import (
    PntSchedulingProblem,
    State,
    Action,
    ServiceWindow,
    UsrId,
    ReqId,
    Request,
)
from .Solver import Solver


class DiscreteTimeIpSolver(Solver):
    def __init__(self, problem: PntSchedulingProblem):
        self.problem: PntSchedulingProblem = problem
        self.solution: np.array = None

    def solve(self, s: State, time_step_factor: float) -> list[tuple[State, Action]]:

        req_dict: dict[ReqId, Request] = {req.id: req for req in self.problem.requests}
        req_by_user: dict[UsrId, list[ReqId]] = {}
        for req in req_dict.values():
            if req.usr_id not in req_by_user:
                req_by_user[req.usr_id] = []
            req_by_user[req.usr_id].append(req.id)

        N_sat = self.problem.N_sat
        N_req = len(req_dict)
        Dt = self.problem.t_final / self.problem.CN0_norm.shape[2]
        dt = time_step_factor * Dt
        N_t = int(self.problem.t_final / dt)
        t_curr = self.problem.current_time

        T = [
            ((req.dur - s.req_times[req_id]) / dt if req.ta <= t_curr else 0)
            for req_id, req in req_dict.items()
        ]
        T = np.ceil(T).astype(int)
        transition_times = np.ceil(self.problem.transition_times / dt).astype(int)
        rewards = np.zeros((N_sat, N_req, N_t))
        service_windows = np.zeros((N_sat, N_req, N_t))
        data_gen_requests = np.zeros((N_sat, N_req, N_t))
        energy_gen_requests = np.zeros((N_sat, N_req, N_t))

        # Service windows
        win_list = (
            self.problem.constr_windows
            if self.problem.constr
            else self.problem.service_windows
        )
        for sat_id in range(N_sat):
            for req_id, req in enumerate(req_dict.values()):
                if req.ta > t_curr:
                    continue
                windows = [
                    win
                    for win in win_list
                    if win.usr_id == req.usr_id and win.sat_id == sat_id
                ]
                starts = np.ceil([w.ts / dt for w in windows]).astype(int)
                ends = np.floor([w.te / dt for w in windows]).astype(int)
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
                    service_windows[sat_id, req_id, ts:te] = 1
                    payload_on = w.usr_id >= 0
                    f = time_step_factor
                    if payload_on:
                        r = self.problem.CN0_norm[
                            sat_id, w.usr_id, (f * ts) : (f * te) : f
                        ]
                        rewards[sat_id, req_id, ts:te] = r * dt
                        assert not np.isnan(rewards).any()
                        data_gen_requests[sat_id, req_id, ts:te] = (
                            self.problem.payload_data_gen * dt
                        )
                        energy_gen_requests[sat_id, req_id, ts:te] = (
                            self.problem.payload_energy_gen * dt
                        )

        energy_gen = np.zeros((N_sat, N_t))
        data_gen = np.zeros((N_sat, N_t))
        for sat_id in range(N_sat):
            for t in range(N_t):
                args = (sat_id, t * dt, (t + 1) * dt)
                energy_gen[sat_id, t] = self.problem.energy_gen_func(*args)
                data_gen[sat_id, t] = self.problem.data_gen_func(*args)

        # Variables
        x = [cp.Variable((N_req, N_t), boolean=True) for _ in range(N_sat)]

        # Objective
        lambda_sep = 1e-5 / (N_req * N_t)
        discount = np.power(0.99, np.arange(10))
        inner_sum = [
            cp.sum(
                cp.multiply(x[i], rewards[i] * discount)
                - lambda_sep * cp.sum(cp.abs(cp.diff(x[i], axis=1)))
            )
            for i in range(N_sat)
        ]
        objective = cp.Maximize(cp.sum(inner_sum))

        # Constraints
        constraints = []

        # Current time
        ti = int(np.ceil(s.t[sat_id] / dt))
        if ti > 0:
            for sat_id in range(N_sat):
                constraints.append(x[sat_id][:, :ti] == 0)

        for sat_id in range(N_sat):
            tmp_data = data_gen[sat_id] + cp.sum(
                cp.multiply(data_gen_requests[sat_id], x[sat_id]), axis=0
            )
            tmp_energy = energy_gen[sat_id] + cp.sum(
                cp.multiply(energy_gen_requests[sat_id], x[sat_id]), axis=0
            )
            for t in range(N_t + 1):
                constraints.append(
                    s.D[sat_id] + cp.sum(tmp_data[:t]) <= self.problem.max_data
                )
                constraints.append(
                    s.E[sat_id] + cp.sum(tmp_energy[:t]) >= self.problem.min_energy
                )

        # One action at a time for each satellite
        constraints.extend([cp.sum(x[sat_id], axis=0) <= 1 for sat_id in range(N_sat)])

        # One satellite at a time for each user
        for req_ids in req_by_user.values():
            for t in range(N_t):
                arr = [
                    x[sat_id][req_id, t]
                    for sat_id in range(N_sat)
                    for req_id in req_ids
                ]
                constraints.append(cp.sum(arr) <= 1)

        # Duration
        inner_sum = [cp.sum(x[sat_id], axis=1) for sat_id in range(N_sat)]
        constraints.append(cp.sum(inner_sum) <= T)

        # Service window
        constraints.extend(
            [x[sat_id] <= service_windows[sat_id] for sat_id in range(N_sat)]
        )
        for sat_id in range(N_sat):
            for i in range(N_req):
                for j in range(N_req):
                    if i == j:
                        continue
                    for t in range(N_t - 1):
                        tt = transition_times[sat_id, i, j]
                        # Transition times
                        constraints.append(
                            x[sat_id][j, t : t + tt + 1] <= 1 - x[sat_id][i, t]
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
        for sat_id in range(N_sat):
            self.solution[sat_id] = np.round(x[sat_id].value).astype(int)

        # Policy
        # policy = self.construct_policy(s, dt)
        # return policy

        # def construct_policy(
        #     self, s: State, time_step: float
        # ) -> list[tuple[State, Action]]:

        # Dictionary mapping request and satellite ids to service windows
        service_windows_dict: dict[int, dict[int, list[ServiceWindow]]] = {
            req.id: {
                sat_id: [
                    w
                    for w in self.problem.service_windows
                    if w.usr_id == req.usr_id and w.sat_id == sat_id
                ]
                for sat_id in range(self.problem.N_sat)
            }
            for req in req_dict.values()
        }

        # Get window id given a time interval and a request id
        def get_window(ts: int, te: int, request_id: int, sat_id: int) -> int:
            windows = [
                win
                for win in service_windows_dict[request_id][sat_id]
                if ts >= int(np.ceil(win.ts / dt)) - 1
                and te <= int(np.floor(win.te / dt)) + 1
            ]
            return windows[0]

        # Create actions by iterating over the requests
        actions: list[Action] = []
        for req_id, req in enumerate(req_dict.values()):
            for sat_id in range(self.problem.N_sat):
                start, end = util.get_start_end_indexes(self.solution[sat_id, req_id])
                for ts, te in zip(start, end):
                    a = Action(sat_id=sat_id, ts=ts * dt, dur=(te - ts) * dt, req=req)
                    actions.append(a)
        actions.sort(key=lambda a: a.ts)

        tmp_d = [
            data_gen[sat_id]
            + np.sum(
                np.multiply(data_gen_requests[sat_id], self.solution[sat_id]), axis=0
            )
            for sat_id in range(N_sat)
        ]
        tmp_e = [
            energy_gen[sat_id]
            + np.sum(
                np.multiply(energy_gen_requests[sat_id], self.solution[sat_id]), axis=0
            )
            for sat_id in range(N_sat)
        ]

        # Create policy
        policy = []
        for a in actions:
            # Fix the duration of the action
            policy.append((s, a))
            s = self.problem.transition_function(s, a)

        policy.append((s, None))
        policy = self.problem.clean_policy(policy)
        return policy
