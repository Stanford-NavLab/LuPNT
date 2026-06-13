import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pylupnt as pnt
import os
from pymoo.indicators.hv import Hypervolume


def plot_coverage_dop(
    coverage_phases, dop_phases, prob_phases, x_user, user_used, tspan, dop_thresholds, perctile=50
):
    """
    Plot the coverage and DOP for each phase.

    Parameters:
    coverage_phases : list of np.ndarray [n_user, lent]
        Coverage for each phase.
    dop_phases : list of np.ndarray [n_user, lent]
        DOP for each phase.
    x_user : np.ndarray [n_users, 3]
        User positions in Cartesian coordinates.
    tspan : np.ndarray [lent]
        Time span for propagation.
    dop_thresholds : list of float
        DOP thresholds for each phase.
    """

    # compute latitudes and longitudes of user positions
    n_phase = len(coverage_phases)
    n_user = x_user.shape[0]
    lent = len(tspan)
    fig, axes = plt.subplots(3, n_phase, figsize=(7.5 * n_phase, 15))

    for phase in range(n_phase):

        latitudes = np.degrees(np.arcsin(x_user[user_used[phase], 2] / pnt.R_MOON))  # [deg]
        longitudes = np.degrees(
            np.arctan2(x_user[user_used[phase], 1], x_user[user_used[phase], 0])
        )  # [deg]

        dop_phase = dop_phases[phase]  # (n_user, lent, case)
        n_user = dop_phase.shape[0]
        prob_phase = prob_phases[phase]
        dop_target_ratio = (
            np.sum(np.sum(dop_phase <= dop_thresholds[phase], axis=0), axis=0) / lent / n_user
        )  # (case)
        wsum_dop_target_ratio = np.sum(prob_phase * dop_target_ratio)

        # generate a color map for coverage
        cov_phase = np.sum(coverage_phases[phase][:, :, 0], axis=1) / lent  # average over time
        dop_phase_perc = np.percentile(
            dop_phases[phase][:, :, 0], perctile, axis=1
        )  # average over time
        dop_phase_ratio = (
            np.sum(dop_phases[phase][:, :, 0] <= dop_thresholds[phase], axis=1) / lent
        )  # ratio of DOP below threshold

        print("----------------------------------------------------")
        print("Phase: ", phase + 1)
        print("  Failure probs: ", prob_phase[prob_phase > 0])
        print("  Total failure prob: ", np.sum(prob_phase))
        print("  DOP Ratio under target (per case): ", dop_target_ratio)
        print("  DOP Ratio under target (no weight): {:.2f}%".format(dop_target_ratio[0] * 100))
        print(
            "  DOP Ratio under target (w/ weight): {:.2f}%".format(wsum_dop_target_ratio * 100),
            "  <------ objective!",
        )

        cmap = plt.get_cmap("viridis")
        norm = plt.Normalize(vmin=0, vmax=np.nanmax(cov_phase))
        ax = axes[0, phase]
        ax.set_title(f"4-Fold Coverage Phase {phase + 1}")
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.grid(True)
        sc = ax.scatter(longitudes, latitudes, c=cov_phase, cmap=cmap, norm=norm, s=10)
        fig.colorbar(sc, ax=axes[0, phase], label="Coverage")

        # generate a color map for DOP
        cmap_dop = plt.get_cmap("viridis")
        norm_dop = plt.Normalize(vmin=0, vmax=20)
        ax = axes[1, phase]
        ax.set_title(f"DOP Phase {phase + 1} (95% values)")
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.grid(True)
        sc_dop = ax.scatter(
            longitudes, latitudes, c=dop_phase_perc, cmap=cmap_dop, norm=norm_dop, s=10
        )
        fig.colorbar(sc_dop, ax=axes[1, phase], label="DOP")

        # DOP ratios
        cmap_dop_ratio = plt.get_cmap("viridis")
        norm_dop_ratio = plt.Normalize(vmin=0, vmax=1)
        ax = axes[2, phase]
        ax.set_title(f"DOP Phase {phase + 1} (Ratio smaller than Threshold)")
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.grid(True)
        sc_dop_ratio = ax.scatter(
            longitudes, latitudes, c=dop_phase_ratio, cmap=cmap_dop_ratio, norm=norm_dop_ratio, s=10
        )
        fig.colorbar(sc_dop_ratio, ax=axes[2, phase], label="DOP Ratio")


def plot_ures(ures, x_user, tspan, n_plot=5):
    """
    Plot user range errors.
    """
    n_users, n_sat, lent = ures.shape

    # select random 5 users
    user_indices = np.random.choice(n_users, size=n_plot, replace=False)
    ures = ures[user_indices]

    rows = int(np.ceil(n_plot / 2))
    cols = 2
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 7.5, rows * 5))

    for i in range(n_plot):
        ax = axes.flatten()[i]
        lat = np.degrees(np.arcsin(x_user[user_indices[i], 2] / pnt.R_MOON))
        lon = np.degrees(np.arctan2(x_user[user_indices[i], 1], x_user[user_indices[i], 0]))
        for j in range(n_sat):
            ax.plot(tspan / 3600, ures[i, j, :], label=f"Sat {j+1}")
        ax.set_title(f"User lat: {lat:.2f}, lon: {lon:.2f}")
        ax.set_xlabel("Time [hr]")
        ax.set_ylabel("Range Error [m]")
        ax.grid(True)

    plt.tight_layout()
    plt.show()


def plot_user_positions(
    x_user, lat_range=(-90, 90), marker_size=2, camera_view=None, fignaeme=None
):
    """
    Plot user positions on the Moon with elevation mask.

    Parameters:
    x_user : np.ndarray [n_users, 3]
        Cartesian coordinates of user positions.
    elev_mask_deg : float
        Elevation mask in degrees.
    """
    fig = go.Figure()
    pnt.plot.plot_body(
        fig,
        pnt.MOON,
        size_factor=2,
        alpha=0.5,
    )

    # plot user positions as red dots
    pnt.plot.scatter(fig, 1.01 * x_user, color="red", marker_size=marker_size)

    # plot the points within the latitude range
    lat_mask = np.logical_and(
        np.degrees(np.arcsin(x_user[:, 2] / pnt.R_MOON)) >= lat_range[0],
        np.degrees(np.arcsin(x_user[:, 2] / pnt.R_MOON)) <= lat_range[1],
    )

    pnt.plot.scatter(fig, 1.01 * x_user[lat_mask], color="blue", marker_size=marker_size)

    # set the fontsize of the axes
    fs = 12
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title="X [km]",
                title_font=dict(size=fs),
                tickfont=dict(size=fs),
            ),
            yaxis=dict(
                title="Y [km]",
                title_font=dict(size=fs),
                tickfont=dict(size=fs),
            ),
            zaxis=dict(
                title="Z [km]",
                title_font=dict(size=fs),
                tickfont=dict(size=fs),
            ),
        ),
    )

    # set the camera view
    if camera_view is not None:
        pnt.plot.set_view(fig, *camera_view)
    else:
        pnt.plot.set_view(fig, -80, 20, 2.5)

    if fignaeme is not None:
        fig.write_image(fignaeme)

    fig.show()


def plot_gs(llas, labels, markersize=50, fontsize=9, figname=None):
    """
    Plot ground station locations on Earth.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Ground Station Locations")

    # load earth topo image
    img_dir = os.path.join(pnt.get_basepath(), "topo", "earth_surface.jpg")
    img = plt.imread(img_dir)
    # flip the image vertically and horizontally
    img = img[::-1, ::-1]
    ax.imshow(img, extent=[-180, 180, -90, 90], alpha=0.2)

    n_gs = llas.shape[0]
    for i in range(n_gs):
        lat, lon, alt = llas[i]
        if lon > 180:
            lon -= 360
        ax.scatter(lon, lat, s=markersize, label=labels[i])
        ax.text(lon + 2, lat + 2, labels[i], fontsize=fontsize)

    if figname is not None:
        plt.savefig(figname)
    plt.show()


def plot_survive_rates(fail_model, launch_years, eval_years, age_max, figname=None):

    n_phase = len(launch_years)

    # ---- Plotting survival function ----
    for i in range(n_phase):
        print("")
        survive_rates = []
        for j in range(i + 1):
            survive_rates.append(fail_model.survival(eval_years[i] - launch_years[j]))
        print("phase: {}".format(i))
        print("  launch year:", launch_years[i])
        print("  eval year:", eval_years[i])
        print("  survive rates: ", survive_rates)

    t = np.linspace(0, age_max, 100)
    S = fail_model.survival(t)
    plt.figure(figsize=(8, 4))
    plt.plot(t, S, label="Survival Function", color="blue")
    plt.title("NavSat Survival Function")
    plt.xlabel("Age (years)")
    plt.ylabel("Survival Probability")
    plt.grid()
    # plt.axhline(0.5, color='red', linestyle='--', label='50% Survival')
    plt.axvline(1, color="purple", linestyle="--", label="Stage 1")
    plt.axvline(6, color="green", linestyle="--", label="Stage 2")
    plt.axvline(11, color="orange", linestyle="--", label="Stage 3")
    plt.ylim(0.5, 1)
    plt.xlim(0, age_max)
    plt.legend()
    if figname is not None:
        plt.savefig(figname)
    plt.show()


def plot_history_hv(res, figname=None):

    n_obj = res.F.shape[1]

    hist = res.history
    n_evals = []
    hist_F = []  # the objective space values in each generation
    hist_cv = []  # constraint violation in each generation
    hist_cv_avg = []  # average constraint violation in the whole population

    for algo in hist:
        # retrieve the optimum from the algorithm
        opt = algo.opt

        # store the number of function evaluations
        n_evals.append(algo.evaluator.n_eval)

        # store the least contraint violation and the average in each population
        hist_cv.append(opt.get("CV").min())
        hist_cv_avg.append(algo.pop.get("CV").mean())

        # filter out only the feasible and append and objective space values
        feas = np.where(opt.get("feasible"))[0]
        hist_F.append(opt.get("F")[feas])

    k = np.where(np.array(hist_cv) <= 0.0)[0].min()
    print(f"At least one feasible solution in Generation {k} after {n_evals[k]} evaluations.")

    metric = Hypervolume(
        ref_point=np.ones(n_obj),
        norm_ref_point=False,
        zero_to_one=True,
        ideal=np.zeros(n_obj),
        nadir=np.ones(n_obj),
    )

    hv = [metric.do(_F) for _F in hist_F]

    fig = plt.figure(figsize=(7, 5))
    plt.plot(n_evals, hv, color="black", lw=0.7, label="Avg. CV of Pop")
    plt.scatter(n_evals, hv, facecolor="none", edgecolor="black", marker="p")
    plt.title("Convergence")
    plt.xlabel("Function Evaluations")
    plt.ylabel("Hypervolume")
    plt.grid(True)
    if figname is not None:
        plt.savefig(figname)
    plt.show()


def plot_coverage_map(
    data_dir, walker_patterns, smas, incs, plot_pole=False, hybrid=False, n_cols=4, with_text=False
):
    # convert smas to km
    smas = smas / 1e3

    ## plot coverage for each walker pattern
    if plot_pole:
        validx = 1
    else:
        validx = 0  # index for coverage metric to plot (0: 4-fold coverage, 1: pole coverage)

    idx_plot = []
    for k, pattern in enumerate(walker_patterns):
        n_sat = pattern[0] * pattern[1]
        if hybrid:
            n_sat *= 2
        idx_plot.append(k)
        # if plot_pole:
        #     if n_sat <= 15:
        #         idx_plot.append(k)
        # else:
        #     if n_sat >= 16:
        #         idx_plot.append(k)

    n_walker = len(idx_plot)
    n_cols = int(min(n_walker, n_cols))
    n_rows = int(np.ceil(n_walker / n_cols))
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))

    k = 0
    for ip, pattern in enumerate(walker_patterns):
        # load data
        datafile = os.path.join(data_dir, f"walker_{pattern[0]}_{pattern[1]}.npy")
        if os.path.exists(datafile):
            results = np.load(datafile)
        else:
            print(f"Data file {datafile} not found. Skipping.")
            continue

        if ip in idx_plot:
            ax = axs[k // n_cols, k % n_cols]
            incd = np.rad2deg(incs)
            igrid, agrid = np.meshgrid(incd, smas)
            coverage = results[:, :, validx] * 100  # convert to percentage
            # contourf plot (0-100)
            cp = ax.contourf(
                incd, smas, coverage, cmap="viridis", levels=50, alpha=0.3, vmin=0, vmax=100
            )
            mp = ax.scatter(igrid, agrid, c=coverage, cmap="viridis", vmin=0, vmax=100)

            cov_over_99 = coverage > 99
            ax.scatter(
                igrid[cov_over_99],
                agrid[cov_over_99],
                c=coverage[cov_over_99],
                cmap="viridis",
                vmin=0,
                vmax=100,
                edgecolor="red",
                s=50,
                label="over 99% coverage",
            )

            # plot the sma limit
            ecc = np.sqrt(1 - 5 / 3 * np.cos(incs) ** 2)
            sma_limit = pnt.R_MOON / (1 - ecc) + 100
            ax.plot(np.rad2deg(incs), sma_limit, "k--")

            if hybrid:
                ax.set_title(
                    f"({pattern[0]} x 2) plane / {pattern[1]} sat (total {pattern[0]*pattern[1]*2} sat)",
                    fontsize=16,
                )
            else:
                ax.set_title(
                    f"{pattern[0]} plane x {pattern[1]} sat (total {pattern[0]*pattern[1]} sat)",
                    fontsize=16,
                )
            ax.set_xlabel("Inclination [deg]", fontsize=16)
            ax.set_ylabel("Semi-major axis [km]", fontsize=16)
            ax.tick_params(axis="both", which="major", labelsize=14)
            ax.set_yticks(smas[::2])
            ax.set_xticks(incd[::2])
            ax.set_ylim(pnt.R_MOON * 1e-3, smas[-1] + 500)
            ax.set_xlim(incd[0] - 0.5, 70.5)
            ax.grid(True)
            # print the text of coverage
            if with_text:
                for i in range(len(smas)):
                    for j in range(len(incd)):
                        eps_h = 200
                        eps_i = 0.5
                        if coverage[i, j] > 0:
                            ax.text(
                                incd[j] + eps_i,
                                smas[i] + eps_h,
                                f"{coverage[i, j]:.2f}",
                                ha="center",
                                va="center",
                                fontsize=9,
                            )
            k += 1
            # colorbar for the counterf plot
            cbar = fig.colorbar(mp, ax=ax)
            cbar.set_label("Coverage (%)", fontsize=14)
            ax.legend(loc="lower right", fontsize=14)

    # tight layout
    fig.tight_layout()

    return fig


def plot_dop_gridsearch_polar(
    sma_polar, incs, walker_patterns_polar, results_polar, is_fault=False
):

    fig, axes = plt.subplots(2, 1, figsize=(8, 6))

    # print the 99% coverage and PDOP ratio satisfcatction
    all_walker_patterns = results_polar.keys()
    min_cov99_sats = np.inf
    min_dop6_sats = np.inf

    for walker_pattern in all_walker_patterns:
        nsat = int(walker_pattern[0] * walker_pattern[1])
        print("Walker pattern: plane{0} x sat{1}".format(walker_pattern[0], walker_pattern[1]))
        if np.max(results_polar[walker_pattern][:, 0]) > 0.99:
            print(
                "reached 99% south pole coverage: {0:.2f}%".format(
                    100 * np.max(results_polar[walker_pattern][:, 0])
                )
            )
            if nsat < min_cov99_sats:
                min_cov99_sats = nsat
        if np.max(results_polar[walker_pattern][:, 2]) > 0.99:
            print(
                "reach 99% south pole PDOP under 6: {0:.2f}%".format(
                    100 * np.max(results_polar[walker_pattern][:, 2])
                )
            )
            if nsat < min_dop6_sats:
                min_dop6_sats = nsat
        if np.max(results_polar[walker_pattern][:, 3]) > 0.9:
            print(
                "reach 90% south pole PDOP under 3: {0:.2f}%".format(
                    100 * np.max(results_polar[walker_pattern][:, 3])
                )
            )

    print("Minimum number of satellites to reach 99% coverage: ", min_cov99_sats)
    print("Minimum number of satellites to reach 99% PDOP under 6: ", min_dop6_sats)

    for k, walker_pattern in enumerate(walker_patterns_polar):
        nsat = int(walker_pattern[0] * walker_pattern[1])

        axes[0].plot(
            np.rad2deg(incs),
            100 * results_polar[walker_pattern][:, 0],
            "o-",
            label="{0} planes x {1} sats = {2} sats".format(
                walker_pattern[0], walker_pattern[1], nsat
            ),
        )
        axes[1].plot(
            np.rad2deg(incs),
            100 * results_polar[walker_pattern][:, 2],
            "o-",
            label="{0} planes x {1} sats = {2} sats".format(
                walker_pattern[0], walker_pattern[1], nsat
            ),
        )

        axes[0].set_xlabel("Inclination [deg]", fontsize=14)
        axes[1].set_xlabel("Inclination [deg]", fontsize=14)
        axes[0].set_ylabel("Pole Coverage (4-fold) [%]", fontsize=14)
        axes[1].set_ylabel("Pole GDOP under 6 ratio [%]", fontsize=14)
        axes[0].grid(True)
        axes[1].grid(True)
        axes[0].legend(fontsize=10, ncol=2)
        axes[1].legend(fontsize=10, ncol=2)
        axes[0].set_ylim(90, 101)
        axes[1].set_ylim(0, 101)
        axes[0].axhline(99, color="k", linestyle="--", alpha=0.5)

    fig.tight_layout()
    filename = "figs/coverage/dop_gridsearch_polar_{0:.0f}.pdf".format(sma_polar * 1e-3)

    if is_fault:
        filename = "figs/coverage/dop_gridsearch_polar_fault_{0:.0f}.pdf".format(sma_polar * 1e-3)

    plt.savefig(filename, dpi=300)
    fig.show()


def plot_dop_gridsearch_hybrid(
    sma_hybrid, incs, walker_patterns_hybrid, results_hybrid, is_fault=False
):

    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    results_idx = [0, 5, 2, 6]
    results_labels = [
        "Global Coverage (4-fold)",
        "South Pole Coverage (4-fold)",
        "Global GDOP under 6 ratio [%]",
        "South Pole GDOP under 6 ratio [%]",
    ]
    axes = axes.flatten()

    all_walker_patterns = results_hybrid.keys()

    min_global_cov99_sats = np.inf
    min_global_dop6_sats = np.inf
    min_pole_cov99_sats = np.inf
    min_pole_dop6_sats = np.inf

    for walker_pattern in all_walker_patterns:
        nsat = int(2 * walker_pattern[0] * walker_pattern[1])
        print(
            "Walker pattern: plane (2 x {0}) x sat{1}".format(walker_pattern[0], walker_pattern[1])
        )
        if np.max(results_hybrid[walker_pattern][:, 0]) > 0.99:
            print(
                "reached 99% global coverage: {0:.2f}%".format(
                    100 * np.max(results_hybrid[walker_pattern][:, 0])
                )
            )
            if nsat < min_global_cov99_sats:
                min_global_cov99_sats = nsat
        if np.max(results_hybrid[walker_pattern][:, 2]) > 0.99:
            print(
                "reach 99% global GDOP under 6: {0:.2f}%".format(
                    100 * np.max(results_hybrid[walker_pattern][:, 2])
                )
            )
            if nsat < min_global_dop6_sats:
                min_global_dop6_sats = nsat
        if np.max(results_hybrid[walker_pattern][:, 5]) > 0.9:
            print(
                "reach 99% south pole coverage: {0:.2f}%".format(
                    100 * np.max(results_hybrid[walker_pattern][:, 5])
                )
            )
            if nsat < min_pole_cov99_sats:
                min_pole_cov99_sats = nsat
        if np.max(results_hybrid[walker_pattern][:, 6]) > 0.99:
            print(
                "reach 99% south pole GDOP under 6: {0:.2f}%".format(
                    100 * np.max(results_hybrid[walker_pattern][:, 6])
                )
            )
            if nsat < min_pole_dop6_sats:
                min_pole_dop6_sats = nsat

    print("Minimum number of satellites to reach 99% global coverage: ", min_global_cov99_sats)
    print("Minimum number of satellites to reach 99% global GDOP under 6: ", min_global_dop6_sats)
    print("Minimum number of satellites to reach 99% south pole coverage: ", min_pole_cov99_sats)
    print("Minimum number of satellites to reach 99% south pole GDOP under 6: ", min_pole_dop6_sats)

    for k, walker_pattern in enumerate(walker_patterns_hybrid):
        nsat = int(2 * walker_pattern[0] * walker_pattern[1])

        for idx, label, row in zip(results_idx, results_labels, range(len(results_idx))):
            axes[row].plot(
                np.rad2deg(incs),
                100 * results_hybrid[walker_pattern][:, idx],
                "o-",
                label="({0} x 2) planes x {1} sats = {2} sats".format(
                    walker_pattern[0], walker_pattern[1], nsat
                ),
            )

    for row in range(len(results_idx)):
        label = results_labels[row]
        axes[row].set_xlabel("Inclination [deg]", fontsize=14)
        axes[row].set_ylabel(label, fontsize=14)
        axes[row].grid(True)
        axes[row].legend(fontsize=10.5, ncol=2)
        # increase font size of axes
        axes[row].tick_params(axis="both", which="major", labelsize=12)
        axes[row].set_xlim(39, 66)
        if "coverage" in label.lower():
            axes[row].set_ylim(90, 101)
            axes[row].axhline(99, color="k", linestyle="--", alpha=0.5)
        else:  # DOP plots
            axes[row].set_ylim(0, 101)

    fig.tight_layout()
    filename = "figs/coverage/dop_gridsearch_hybrid_{0:.0f}.pdf".format(sma_hybrid * 1e-3)

    if is_fault:
        filename = "figs/coverage/dop_gridsearch_hybrid_fault_{0:.0f}.pdf".format(sma_hybrid * 1e-3)

    plt.savefig(filename, dpi=300)
    fig.show()


def plot_dop_gridsearch_circular(
    sma_circular, incs, walker_patterns_circular, results_circular, is_fault=False
):

    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    results_idx = [0, 5, 2, 6]
    results_labels = [
        "Global Coverage (4-fold)",
        "South Pole Coverage (4-fold)",
        "Global GDOP under 6 ratio [%]",
        "South Pole GDOP under 6 ratio [%]",
    ]
    axes = axes.flatten()

    all_walker_patterns = results_circular.keys()
    min_global_cov99_sats = float("inf")
    min_global_dop6_sats = float("inf")
    min_pole_cov99_sats = float("inf")
    min_pole_dop6_sats = float("inf")

    sma_24hr = 9750.7 * 1000  # km, the sma for 24hr period
    sma_24hr_index = int(np.where(sma_circular < sma_24hr)[0][-1])
    print("SMA: ", sma_circular[sma_24hr_index] * 1e-3)

    for walker_pattern in all_walker_patterns:
        nsat = int(walker_pattern[0] * walker_pattern[1])
        print("Walker pattern: plane {0} x sat {1}".format(walker_pattern[0], walker_pattern[1]))
        if np.max(results_circular[walker_pattern][:sma_24hr_index, 0]) > 0.99:
            print(
                "reached 99% global coverage: {0:.2f}%".format(
                    100 * np.max(results_circular[walker_pattern][:sma_24hr_index, 0])
                )
            )
            if nsat < min_global_cov99_sats:
                min_global_cov99_sats = nsat
        if np.max(results_circular[walker_pattern][:sma_24hr_index, 2]) > 0.99:
            print(
                "reach 99% global GDOP under 6: {0:.2f}%".format(
                    100 * np.max(results_circular[walker_pattern][:sma_24hr_index, 2])
                )
            )
            if nsat < min_global_dop6_sats:
                min_global_dop6_sats = nsat
        if np.max(results_circular[walker_pattern][:sma_24hr_index, 5]) > 0.9:
            print(
                "reach 90% south pole coverage: {0:.2f}%".format(
                    100 * np.max(results_circular[walker_pattern][:sma_24hr_index, 5])
                )
            )
            if nsat < min_pole_cov99_sats:
                min_pole_cov99_sats = nsat
        if np.max(results_circular[walker_pattern][:sma_24hr_index, 6]) > 0.99:
            print(
                "reach 99% south pole GDOP under 6: {0:.2f}%".format(
                    100 * np.max(results_circular[walker_pattern][:sma_24hr_index, 6])
                )
            )
            if nsat < min_pole_dop6_sats:
                min_pole_dop6_sats = nsat

    print("Minimum number of satellites to reach 99% global coverage: ", min_global_cov99_sats)
    print("Minimum number of satellites to reach 99% global GDOP under 6: ", min_global_dop6_sats)
    print("Minimum number of satellites to reach 90% south pole coverage: ", min_pole_cov99_sats)
    print("Minimum number of satellites to reach 99% south pole GDOP under 6: ", min_pole_dop6_sats)

    for k, walker_pattern in enumerate(walker_patterns_circular):
        nsat = int(walker_pattern[0] * walker_pattern[1])

        for idx, label, row in zip(results_idx, results_labels, range(len(results_idx))):
            axes[row].plot(
                sma_circular * 1e-3,
                100 * results_circular[walker_pattern][:, idx],
                "o-",
                label="{0} planes x {1} sats = {2} sats".format(
                    walker_pattern[0], walker_pattern[1], nsat
                ),
            )

        if np.max(results_circular[walker_pattern][:, idx]) > 0.99:
            print(
                "reached 99% {0:s}: {1:.2f}%".format(
                    label, 100 * np.max(results_circular[walker_pattern][:, idx])
                )
            )

    for row in range(len(results_idx)):
        label = results_labels[row]
        axes[row].set_xlabel("Semi-Major Axis [km]", fontsize=14)
        axes[row].set_ylabel(label, fontsize=14)
        axes[row].grid(True)
        axes[row].legend(fontsize=10.5, ncol=2)
        # increase font size of axes
        axes[row].tick_params(axis="both", which="major", labelsize=12)
        # axes[row].set_xlim(39, 66)
        if "coverage" in label.lower():
            axes[row].set_ylim(90, 101)
            axes[row].axhline(99, color="k", linestyle="--", alpha=0.5)
        else:  # DOP plots
            axes[row].set_ylim(0, 101)

    fig.tight_layout()
    if is_fault:
        plt.savefig("figs/coverage/dop_gridsearch_circular_fault.pdf", dpi=300)
    else:
        plt.savefig("figs/coverage/dop_gridsearch_circular.pdf", dpi=300)

    fig.show()
