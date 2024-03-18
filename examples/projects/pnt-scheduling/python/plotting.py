import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
import pylupnt as pnt
from matplotlib.colors import TABLEAU_COLORS
from algorithms import Request, ServiceWindow, State, Action, PntSchedulingProblem
from typing import Tuple


def plot_satellites_users(
    rv_moon_sat, rv_moon_user, rv_moon_earth, rv_moon_sun, user_type
) -> None:
    # Satellite orbit in MI frame
    fig = pnt.plots.Plot3D(figsize=(5, 5), elev=15, azim=-40)
    fig.plot_surface(pnt.MOON, scale=3)
    fig.label_axis()
    rv_surface = rv_moon_user[user_type == "surface", :, :]
    rv_orbital = rv_moon_user[user_type == "orbital", :, :]
    fig.scatter(
        rv_surface[:, 0, :3],
        color="tab:red",
        mask=False,
    )
    fig.plot(
        rv_orbital[:, :, :3],
        color="tab:red",
        mask=True,
    )
    fig.ax.plot([], [], color="tab:red", label="Users")
    fig.plot(rv_moon_sat[:, :3], label="Satellite", color="tab:blue", mask=False)

    # Earth and Sun positions
    earth_dir = rv_moon_earth[0, :3]
    earth_dir = earth_dir / np.linalg.norm(earth_dir) * 3 * pnt.R_MOON
    sun_dir = rv_moon_sun[0, :3]
    sun_dir = sun_dir / np.linalg.norm(sun_dir) * 3 * pnt.R_MOON
    org = np.zeros(3)
    fig.ax.quiver(
        *org, *earth_dir, color="tab:green", arrow_length_ratio=0.1, label="Earth"
    )
    fig.ax.quiver(
        *org, *sun_dir, color="tab:orange", arrow_length_ratio=0.1, label="Sun"
    )

    # Change z-axis label margin
    fig.ax.xaxis.labelpad = 0
    fig.ax.yaxis.labelpad = 0
    fig.ax.zaxis.labelpad = 0
    fig.ax.tick_params(axis="both", which="major", pad=0)

    fig.ax.set_xlim(-5e3, 5e3)
    fig.ax.set_ylim(-5e3, 5e3)
    fig.ax.set_zlim(-10e3, 5e3)
    fig.ax.set_xticklabels([f"{int(x/1000)}" for x in fig.ax.get_xticks()])
    fig.ax.set_yticklabels([f"{int(y/1000)}" for y in fig.ax.get_yticks()])
    fig.ax.set_zticklabels([f"{int(z/1000)}" for z in fig.ax.get_zticks()])
    fig.ax.set_xlabel("X [$10^3$ km]")
    fig.ax.set_ylabel("Y [$10^3$ km]")
    fig.ax.set_zlabel("Z [$10^3$ km]")
    fig.ax.set_aspect("equal")
    plt.legend()


def plot_elevation_range(
    tspan: np.array, az_el_rho: np.array, user_visible: np.array, axs: plt.Axes = None
) -> None:
    if axs is None:
        fig, axs = plt.subplots(2, 1, figsize=(8, 6))

    # Elevation and range over time
    x = tspan / pnt.SECS_PER_HOUR
    x_lim = (0, tspan[-1] / pnt.SECS_PER_HOUR)

    plt.sca(axs[0])
    el = np.rad2deg(az_el_rho[:, :, 1].copy())
    plt.plot(x, el.T, label="Elevation", lw=0.5)
    el[~user_visible] = np.nan
    plt.plot(x, el.T, lw=3)
    plt.xlabel("Time (h)")
    plt.ylabel("Elevation (deg)")
    plt.xlim(x_lim)
    plt.grid()

    plt.sca(axs[1])
    rho = az_el_rho[:, :, 2].copy()
    plt.plot(x, rho.T, label="Range", lw=0.5)
    rho[~user_visible] = np.nan
    plt.plot(x, rho.T, lw=3)
    plt.xlabel("Time (h)")
    plt.ylabel("Range (km)")
    plt.ylim(0, 15e3)
    plt.xlim(x_lim)
    plt.grid()
    plt.tight_layout()


def plot_service_windows(
    contact_start_ends: list[np.array], duration_factor: float, time_multiplier: float
) -> None:
    plt.figure(figsize=(8, 3))
    for i, start_ends in enumerate(contact_start_ends, 1):
        start_ends = start_ends * time_multiplier
        plt.plot(
            start_ends.T,
            i * np.ones_like(start_ends.T),
            "black",
            lw=2,
        )
        mid = (start_ends[:, 0] + start_ends[:, 1]) / 2
        dur = start_ends[:, 1] - start_ends[:, 0]
        per = dur * duration_factor
        for m, p in zip(mid, per):
            plt.fill_between(
                [m - p / 2, m + p / 2], i - 0.3, i + 0.3, alpha=0.7, color="tab:blue"
            )

    plt.plot([], [], color="black", lw=2, label="Visibility windows")
    plt.fill_between(
        [], [], alpha=0.7, color="tab:blue", label="Average duration per window"
    )
    plt.yticks(np.arange(len(contact_start_ends)) + 1)
    plt.xlabel("Time [h]")
    plt.ylabel("User")
    plt.grid()
    plt.legend(loc="upper right", facecolor="white", framealpha=1)
    plt.title("Service Windows")
    plt.tight_layout()


def plot_requests_service_windows(
    requests: list[Request],
    service_windows: list[ServiceWindow],
    policy: list[tuple[State, Action]] = None,
    ax: plt.Axes = None,
) -> None:
    request_dict: dict[int, Request] = {r.id: r for r in requests}
    window_dict: dict[int, ServiceWindow] = {w.id: w for w in service_windows}
    total_contact: dict[int, float] = {r.id: 0 for r in requests}
    for w in service_windows:
        total_contact[w.request_id] += w.end - w.start

    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)

    # Plot service windows
    for w in service_windows:
        y = w.request_id
        plt.plot(
            [w.start, w.end],
            [y, y],
            "black",
            lw=2,
        )

    if policy is None:
        # Plot average duration per window
        for w in service_windows:
            y = w.request_id
            d = (
                request_dict[w.request_id].duration
                * (w.end - w.start)
                / total_contact[w.request_id]
            )
            m = (w.start + w.end) / 2
            plt.fill_between(
                [m - d / 2, m + d / 2],
                y - 0.3,
                y + 0.3,
                alpha=0.7,
                color="tab:blue",
            )
    else:
        # Plot policy
        for s, a in policy[:-1]:
            y = window_dict[a.window.id].request_id
            plt.fill_between(
                [a.start, a.start + a.duration],
                y - 0.3,
                y + 0.3,
                alpha=0.7,
                color="tab:blue",
                edgecolor=None,
            )

    plt.plot([], [], color="black", lw=2, label="Service windows")
    plt.fill_between([], [], alpha=0.7, label="Average duration per window")
    plt.yticks(np.arange(min(request_dict.keys()), max(request_dict.keys()) + 1))
    plt.xlim(0, np.max([w.end for w in service_windows]))
    plt.xlabel("Time")
    plt.ylabel("Request ID")
    plt.grid()
    plt.gca().invert_yaxis()
    plt.tight_layout()


def plot_resources(
    problem: PntSchedulingProblem,
    policy: list[tuple[State, Action]],
    ax: plt.Axes = None,
) -> None:

    times = []
    data = []
    energy = []
    for s, a in policy:
        times.append(s.time)
        data.append(s.data)
        energy.append(s.energy)
        if a is not None:
            times.append(a.start)
            data.append(
                max(s.data + problem.data_gen_func(s.time, a.start), problem.min_data)
            )
            energy.append(
                min(
                    s.energy + problem.energy_gen_func(s.time, a.start),
                    problem.max_energy,
                )
            )

    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(8, 6))

    plt.sca(ax[0])
    plt.plot(times, data, "tab:blue", label="Data", lw=2)
    plt.hlines(problem.min_data, 0, np.max(times), colors="tab:blue", linestyles="--")
    plt.hlines(problem.max_data, 0, np.max(times), colors="tab:blue", linestyles="--")
    plt.xlabel("Time")
    plt.ylabel("Data [MB]")
    plt.xlim(0, np.max(times))
    plt.grid()

    plt.sca(ax[1])
    plt.plot(times, energy, "tab:green", label="Energy", lw=2)
    plt.hlines(
        problem.min_energy, 0, np.max(times), colors="tab:green", linestyles="--"
    )
    plt.hlines(
        problem.max_energy, 0, np.max(times), colors="tab:green", linestyles="--"
    )
    plt.xlabel("Time")
    plt.ylabel("Energy [kWh]")
    plt.xlim(0, np.max(times))
    plt.grid()

    plt.tight_layout()


def plot_windows(
    targets: np.array,  # Target ids
    vtws: np.array,  # Visibility time windows (start, end)
    durations: np.array,  # Task durations
    #
    task_otws: np.array,  # Opportunity time windows (start, end)
    task_idxs: np.array,  # Target ids
    #
    rewards: np.array = None,  # Rewards
    selected_tasks: np.array = None,  # Selected tasks
    plot_labels: bool = True,
):
    color = TABLEAU_COLORS[0]
    if selected_tasks is None:
        selected_targets = None
    else:
        assert selected_tasks.dtype == bool
        assert task_idxs.dtype == int
        selected_targets = task_idxs[selected_tasks]

    # Sort tasks by start time
    indices = sorted(range(len(vtws)), key=lambda i: vtws[i][0])

    # Add tasks to the minimum number of rows for plotting
    rows = []
    for idx in indices:
        for row in rows:
            if vtws[idx][0] > vtws[row[-1]][1]:
                row.append(idx)
                break
        else:
            rows.append([idx])

    # Plot tasks
    for r, row in enumerate(rows):

        def add_text(ts, te, txt):
            if plot_labels:
                plt.text((ts + te) / 2, r + 0.15, txt, ha="center", va="center")

        def add_rect(ts, te):
            eps = 0
            plt.fill_between(
                [ts + eps, te - eps], r - 0.3, r + 0.3, alpha=0.5, color=color
            )

        def add_empty_rect(ts, te):
            plt.fill_between(
                [ts, te],
                r - 0.3,
                r + 0.3,
                alpha=0.15,
                hatch="/",
                edgecolor="black",
                color=color,
            )

        for idx in row:
            # Opportunity window
            vtw_s = vtws[idx][0]
            vtw_e = vtws[idx][1]
            dur = durations[idx]

            if True:
                # Plot lines for visibility time windows
                plt.hlines(r, vtw_s, vtw_e, colors="k")
                plt.vlines(vtw_s, r - 0.15, r + 0.15, colors="k")
                plt.vlines(vtw_e, r - 0.15, r + 0.15, colors="k")

                # Plot lines for opportunity time windows
                otws = task_otws[task_idxs == idx]
                if len(otws) > 0:
                    for otw in otws:
                        plt.vlines(otw[0], r - 0.05, r + 0.05, colors="k")

                # Plot rectangles for tasks
                if selected_targets is None:
                    # Start times not resolved
                    ts = (vtw_s + vtw_e - dur) / 2
                    te = (vtw_s + vtw_e + dur) / 2
                    add_rect(ts, te)
                    add_text(vtw_s, vtw_e, str(targets[idx]))
                else:

                    if idx not in selected_targets:
                        # Start times resolved and task is not selected
                        ts = (vtw_s + vtw_e - dur) / 2
                        te = (vtw_s + vtw_e + dur) / 2
                        add_empty_rect(ts, te)
                        add_text(vtw_s, vtw_e, str(targets[idx]))
                    else:
                        # Start times resolved and task is selected
                        for ts, te in task_otws[selected_tasks & (task_idxs == idx)]:
                            add_rect(ts, te)
                            add_text(ts, te, str(targets[idx]))

            else:
                # Sun pointing or downlink opportunity
                if start_end_times is None:
                    add_rect(vtw_start, t_end)
                    add_text(vtw_start, t_end, str(opp.id))
                elif opp.id in start_end_times:
                    for ts, te in zip(*start_end_times[opp.id]):
                        add_rect(ts, te)
                        add_text(ts, te, str(opp.id))
                        add_text(ts, te, str(opp.id))

    plt.xlabel("Time")
    plt.xlim(0, np.max(vtws[:, 1]))
    # plt.gca().spines[["left", "top", "right"]].set_visible(False)
    # if start_end_times is not None:
    #     legend = [
    #         plt.Line2D(
    #             [0], [0], color=COLORS[tt], lw=5, label=MODE_NAMES[tt], alpha=0.5
    #         )
    #         for tt in TaskType
    #         if tt not in [TaskType.START]
    #     ]
    #     # Legend with 3 columns at the bottom outside the plot
    #     plt.legend(
    #         handles=legend,
    #         loc="upper center",
    #         bbox_to_anchor=(0.5, 1.3),
    #         ncol=3,
    #     )
    plt.yticks([])
    plt.ylabel("Tasks")
