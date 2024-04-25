import numpy as np
import pprint
import matplotlib.pyplot as plt
import plotting as plotting
from solvers import (
    SmdpForwardSearchSolver,
    SmdpMctsSolver,
    DiscreteTimeIpSolver,
    RuleBasedSolver,
)
from problem import (
    State,
    Action,
    Request,
    ServiceWindow,
    PntSchedulingProblem,
    reset_id_counters,
)
from problem import UsrId, SatId, ReqId, Time, Energy, Data

pp = pprint.PrettyPrinter()
figures_path = "../figures/AA229/"


def main():
    reset_id_counters()
    requests = [
        Request(usr_id=0, ts=0.0, te=10.0, dur=2.0, ta=0.0),
        Request(usr_id=1, ts=3.0, te=10.0, dur=3.0, ta=0.0),
        Request(usr_id=1, ts=0.0, te=10.0, dur=2.0, ta=3.0),
    ]
    request_dict = {r.id: r for r in requests}

    service_windows = [
        # Satellite 0
        ServiceWindow(usr_id=0, sat_id=0, ts=0.0, te=5.0),
        ServiceWindow(usr_id=1, sat_id=0, ts=0.0, te=3.0),
        ServiceWindow(usr_id=1, sat_id=0, ts=6.0, te=10.0),
        # Satellite 1
        ServiceWindow(usr_id=0, sat_id=1, ts=3.0, te=7.0),
        ServiceWindow(usr_id=1, sat_id=1, ts=5.0, te=7.0),
    ]

    N_sat = 2
    N_usr = 2
    N_req = len(requests)
    N_win = len(service_windows)

    transition_times = np.ones((len(service_windows), len(service_windows)))
    transition_times = np.ones((N_sat, N_req, N_req))
    for i in range(N_sat):
        transition_times[i, np.diag_indices(N_req)] = 0
    transition_times[0, :] = 0
    transition_times[:, 0] = 0

    t_step = 1
    t_final = 10
    N_steps = int(t_final / t_step)
    CN0 = np.ones((2, N_usr, N_steps))

    data_gen = -0.5 * np.ones((N_sat, N_steps))
    energy_gen = 0.5 * np.ones((N_sat, N_steps))

    problem = PntSchedulingProblem(
        t_step=t_step,
        t_final=t_final,
        requests=requests,
        service_windows=service_windows,
        transition_times=transition_times,
        CN0=CN0,
        max_energy=2 / 0.2,
        min_energy=2,
        max_data=8,
        min_data=8 / 0.8 * 0.2,
        payload_data_gen=1,
        payload_energy_gen=-1,
        energy_gen=energy_gen,
        data_gen=data_gen,
    )

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))
    plotting.plot_requests_service_windows(requests, service_windows, ax=ax)
    plt.tight_layout()

    # *************************************************************************

    # Parameters
    d = 5
    gamma = 0.7
    N_max = 5
    d_min = 2

    # Config
    problem.set_current_time(0)
    problem.set_current_policy(None)

    # Solve
    s = problem.initial_state()
    solver = SmdpForwardSearchSolver(problem)
    policy = solver.solve(s, d=d, gamma=gamma, N_max=N_max, d_min=d_min)

    # Metrics
    print(f"Total reward: {problem.total_reward(policy, gamma=gamma):.2f}")
    print("Percentage of requests served:")
    percentage = problem.percentage_completed(policy)
    pp.pprint(percentage)
    print("Total:", round(problem.total_reward(policy, gamma=gamma), 2))

    actions = [
        Action(sat_id=0, request=requests[1], ts=6.0, dur=1.0),
        Action(sat_id=1, request=requests[0], ts=4.0, dur=1.0),
    ]

    s = problem.initial_state()
    policy = list[tuple[State, Action]]()
    for a in actions:
        policy.append((s, a))
        s = problem.transition_function(s, a)
    policy.append((s, None))

    fig, axs = plt.subplots(3, 1, figsize=(6, 6))
    print(f"Total reward: {problem.total_reward(policy, gamma=gamma):.2f}")
    plotting.plot_requests_service_windows(requests, service_windows, policy, ax=axs[0])
    plotting.plot_resources(problem, policy, ax=axs[1:])

    # *************************************************************************
    # First pass

    # Config
    old_policy = policy
    problem.set_current_time(3)
    problem.set_current_policy(old_policy, constr=True)
    problem.initial_state()
    fig, axs = plt.subplots(3, 1, figsize=(6, 6))
    print(f"Total reward: {problem.total_reward(policy, gamma=gamma):.2f}")
    plotting.plot_requests_service_windows(
        requests, problem.constr_windows, policy, ax=axs[0]
    )
    plotting.plot_resources(problem, policy, ax=axs[1:])

    # Forward search
    d = 6
    gamma = 0.7
    N_max = 5
    d_min = 10
    problem.set_current_time(3)
    problem.set_current_policy(old_policy, constr=True)
    s = problem.initial_state()

    solver = SmdpForwardSearchSolver(problem)
    policy = solver.solve(s, d=d, gamma=gamma, N_max=N_max, d_min=d_min)
    # Metrics
    print(f"Total reward: {problem.total_reward(policy, gamma=gamma):.2f}")
    print("Percentage of requests served:")
    percentage = problem.percentage_completed(policy)
    pp.pprint(percentage)
    print("Total:", round(problem.total_reward(policy, gamma=gamma), 2))
    fig, axs = plt.subplots(3, 1, figsize=(6, 6))
    print(f"Total reward: {problem.total_reward(policy, gamma=gamma):.2f}")
    plotting.plot_requests_service_windows(
        requests, service_windows, policy, ax=axs[0], current_time=3
    )
    plotting.plot_resources(problem, policy, ax=axs[1:])

    # *************************************************************************

    merged_policy = problem.merge_policies(policy)

    fig, axs = plt.subplots(3, 1, figsize=(6, 6))
    print(f"Total reward: {problem.total_reward(merged_policy, gamma=gamma):.2f}")
    plotting.plot_requests_service_windows(
        requests, service_windows, policy, ax=axs[0], current_time=problem.current_time
    )
    plotting.plot_resources(problem, merged_policy, ax=axs[1:])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
    plt.show()
    print("Done!")
