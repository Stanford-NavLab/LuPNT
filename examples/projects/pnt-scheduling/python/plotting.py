import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
import pylupnt as pnt
from matplotlib.colors import TABLEAU_COLORS
from problem import (
    Request,
    ServiceWindow,
    State,
    Action,
    PntSchedulingProblem,
)
from typing import Tuple

COLORS = [k for k in TABLEAU_COLORS.keys()]
START_IDX = 0


def plot_satellites_users(
    rv_moon_sat: np.ndarray,
    rv_moon_user: np.ndarray,
    rv_moon_earth: np.ndarray,
    rv_moon_sun: np.ndarray,
    user_type: np.ndarray,
) -> None:
    # Satellite orbit in MI frame
    fig = pnt.plots.Plot3D(figsize=(5, 5), elev=15, azim=-50)
    fig.plot_surface(pnt.MOON, scale=3)
    rv_surface = rv_moon_user[user_type == "surface", :, :]
    rv_orbital = rv_moon_user[user_type == "orbital", :, :]
    fig.scatter(
        rv_surface[:, 0, :3],
        color="tab:red",
        mask=False,
        depthshade=False,
        # marker="o",
        s=25,
    )
    fig.plot(
        rv_orbital[:, :, :3],
        color="tab:red",
        mask=True,
        lw=1,
    )
    fig.scatter(
        rv_orbital[:, 0, :3],
        color="tab:red",
        mask=True,
        depthshade=False,
        s=25,
    )
    fig.ax.plot(
        [],
        [],
        color="tab:red",
        label="Users",
        marker="o",
        lw=1,
    )
    fig.ax.plot([], [], color="black", label="Satellites", marker="o", lw=1)

    if len(rv_moon_sat.shape) == 2:
        # Single satellite
        fig.plot(rv_moon_sat[:, :3], color="black", mask=False, lw=1)
        fig.scatter(rv_moon_sat[0, :3], color="black", mask=False, depthshade=False)
    else:
        for i in range(rv_moon_sat.shape[0]):
            fig.plot(rv_moon_sat[i, :, :3], color="black", mask=False, lw=1)
        fig.scatter(
            rv_moon_sat[:, 0, :3], color="black", mask=False, s=25, depthshade=False
        )

    # Earth and Sun positions
    earth_dir = rv_moon_earth[0, :3]
    earth_dir = earth_dir / np.linalg.norm(earth_dir) * 3 * pnt.R_MOON
    sun_dir = rv_moon_sun[0, :3]
    sun_dir = sun_dir / np.linalg.norm(sun_dir) * 3 * pnt.R_MOON
    org = np.zeros(3)
    # fig.ax.quiver(
    #     *org, *earth_dir, color="tab:green", arrow_length_ratio=0.1, label="Earth"
    # )
    # fig.ax.quiver(
    #     *org, *sun_dir, color="tab:orange", arrow_length_ratio=0.1, label="Sun"
    # )

    fig.set_tickpad(0)
    fig.set_tick_multiplier(1e-3)
    fig.set_labels("X [$10^3$ km]", "Y [$10^3$ km]", "Z [$10^3$ km]")
    fig.set_pane_color([1, 1, 1])
    fig.set_labelpad(0, 0, 0)
    fig.set_lims([-5e3, 5e3], [-5e3, 5e3], [-10e3, 5e3])
    plt.legend(facecolor="white", framealpha=1, loc="upper right")


def plot_elevation_range(
    tspan: np.array, az_el_rho: np.array, user_visible: np.array, axs: plt.Axes = None
) -> None:
    if axs is None:
        fig, axs = plt.subplots(2, 1, figsize=(8, 6))

    # Elevation and range over time
    x = tspan / pnt.SECS_PER_HOUR
    x_lim = (0, tspan[-1] / pnt.SECS_PER_HOUR)

    color = "tab:blue"
    plt.sca(axs[0])
    el = np.rad2deg(az_el_rho[:, :, 1].copy())
    plt.plot(x, el.T, lw=0.5, color=color, linestyle="--")
    el[~user_visible] = np.nan
    plt.plot(x, el.T, lw=1.5, color=color)
    plt.xlabel("Time [h]")
    plt.ylabel("Elevation [deg]")
    plt.xlim(x_lim)
    plt.grid()
    plt.plot([], [], lw=1.5, color=color, label="User visible")
    plt.plot([], [], lw=0.5, color=color, linestyle="--", label="User not visible")
    plt.legend(loc="lower right", facecolor="white", framealpha=1)

    plt.sca(axs[1])
    rho = az_el_rho[:, :, 2].copy() / 1e3
    plt.plot(x, rho.T, lw=0.5, color=color, linestyle="--")
    rho[~user_visible] = np.nan
    plt.plot(x, rho.T, lw=1.5, color=color)
    plt.xlabel("Time [h]")
    plt.ylabel("Range [$10^3$ km]")
    plt.ylim(0, 15e3)
    plt.xlim(x_lim)
    plt.grid()
    plt.plot([], [], lw=1.5, color=color, label="User visible")
    plt.plot([], [], lw=0.5, color=color, linestyle="--", label="User not visible")
    plt.legend(loc="lower right", facecolor="white", framealpha=1)

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
    request_dict: dict[int, Request] = {req.id: req for req in requests}
    total_contact: dict[int, float] = {req.id: 0 for req in requests}
    for win in service_windows:
        total_contact[win.user_id] += win.end - win.start

    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)

    # Plot service windows
    N_satellites = len(set([w.satellite_id for w in service_windows]))
    dy = 0.06
    for win in service_windows:
        y = win.user_id + START_IDX
        y += -dy / 2 * N_satellites + dy * N_satellites * win.satellite_id / (
            N_satellites - 1
        )
        plt.plot(
            [win.start, win.end],
            [y, y],
            COLORS[win.satellite_id],
            lw=3,
        )

    dx = 0.01

    if policy is None:
        # Plot average duration per window
        for win in service_windows:
            y = win.user_id + START_IDX
            # y += -0.15 + 0.3 * w.satellite_id / (N_satellites - 1)
            d = (
                request_dict[win.user_id].duration
                * (win.end - win.start)
                / total_contact[win.user_id]
                / N_satellites
            )
            m = (win.start + win.end) / 2
            plt.fill_between(
                [m - d / 2 + dx, m + d / 2 - dx],
                y - 0.4,
                y + 0.4,
                alpha=0.5,
                color=COLORS[win.satellite_id],
            )
    else:
        # Plot policy
        for s, a in policy[:-1]:
            if a.request is None:
                continue
            y = a.request.user_id + START_IDX
            plt.fill_between(
                [a.start + dx, a.start + a.duration - dx],
                y - 0.3,
                y + 0.3,
                alpha=0.5,
                color=COLORS[a.satellite_id],
                edgecolor=None,
            )
            plt.text(
                a.start + a.duration / 2,
                y,
                f"{a.request.id + START_IDX}",
                ha="center",
                va="center",
                color="white",
            )

    plt.plot([], [], color="black", lw=3, label=f"Window")
    plt.plot([], [], color="black", lw=10, label=f"Service", alpha=0.5)
    for i in range(N_satellites):
        plt.plot([], [], color=COLORS[i], lw=3, label=f"Satellite {i+START_IDX}")

    # if policy is None:
    # plt.fill_between([], [], alpha=0.7, label="Average duration")
    plt.yticks(
        np.arange(
            min(request_dict.keys()) + START_IDX,
            max(request_dict.keys()) + 0.5 + START_IDX,
        )
    )
    plt.xlim(0, np.max([w.end for w in service_windows]))
    plt.xlabel("Time")
    plt.ylabel("User")
    plt.ylim(-0.5 + START_IDX, max(request_dict.keys()) - 0.5 + START_IDX)
    plt.grid()
    # legend outside top
    plt.legend(
        facecolor="white",
        framealpha=1,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.3),
        ncol=4,
        frameon=False,
        handlelength=1,
    )
    plt.gca().invert_yaxis()
    plt.tight_layout()


def plot_resources(
    problem: PntSchedulingProblem,
    policy: list[tuple[State, Action]],
    ax: plt.Axes = None,
) -> None:
    times = [[] for _ in range(problem.N_satellites)]
    data = [[] for _ in range(problem.N_satellites)]
    energy = [[] for _ in range(problem.N_satellites)]
    for s, a in policy:
        if a is None:
            continue
        i = a.satellite_id
        times[i].append(s.times[i])
        data[i].append(s.data[i])
        energy[i].append(s.energy[i])
        assert not np.isnan(s.energy).all()
        d = max(
            s.data[i] + problem.data_gen_func(s.times[i], a.start), problem.min_data
        )
        e = min(
            s.energy[i] + problem.energy_gen_func(s.times[i], a.start),
            problem.max_energy,
        )
        assert not np.isnan(e).all()
        times[i].append(a.start)
        data[i].append(d)
        energy[i].append(e)

    for i in range(problem.N_satellites):
        tt = np.linspace(s.times[i], problem.tf, 100)
        d = np.maximum(
            problem.min_data, s.data[i] + problem.data_gen_func(s.times[i], tt)
        )
        e = np.minimum(
            problem.max_energy,
            s.energy[i] + problem.energy_gen_func(s.times[i], tt),
        )
        assert not np.isnan(e).all()
        times[i].extend(tt)
        data[i].extend(d)
        energy[i].extend(e)

    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(8, 6))

    d_lim = problem.max_data / 0.8
    plt.sca(ax[0])
    for i in range(problem.N_satellites):
        x = np.array(times[i])
        y = np.array(data[i]) / d_lim * 100
        plt.plot(x, y, lw=2, color=COLORS[i], label=f"Satellite {i+START_IDX}")
    y = 100 * 0.8
    # plt.hlines(problem.min_data, 0, problem.tf, colors="tab:blue", linestyles="--")
    plt.hlines(
        y, 0, problem.tf, colors="black", linestyles="--", label="Resource constraint"
    )
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Data [\\%]")
    plt.xlim(0, problem.tf)
    plt.ylim(0, 100)
    plt.legend(
        loc="upper center",
        facecolor="white",
        framealpha=1,
        ncol=3,
        bbox_to_anchor=(0.5, 1.4),
        frameon=False,
        handlelength=1,
    )
    plt.grid()

    e_lim = problem.min_energy / 0.2
    plt.sca(ax[1])
    for i in range(problem.N_satellites):
        x = np.array(times[i])
        y = np.array(energy[i]) / e_lim * 100
        plt.plot(x, y, lw=2, color=COLORS[i])
    y = problem.min_energy / e_lim * 100
    plt.hlines(y, 0, problem.tf, colors="black", linestyles="--", label="Min. energy")
    # plt.hlines(problem.max_energy, 0, problem.tf, colors="tab:green", linestyles="--")
    plt.xlabel("Time")
    plt.ylabel("Energy [\\%]")
    plt.xlim(0, problem.tf)
    plt.ylim(0, 100)
    # plt.legend(loc="upper right", facecolor="white", framealpha=1)
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
