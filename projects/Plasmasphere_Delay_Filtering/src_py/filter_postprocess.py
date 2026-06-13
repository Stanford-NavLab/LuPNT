import numpy as np
import matplotlib.pyplot as plt
import h5py


def compute_sise_pos_vel(h5file, eval_tidx, is_smoother=False):

    true_rva_vals = h5file["rva_true_sat"]
    true_clock_vals = h5file["clk_true_sat"]

    if is_smoother:
        est_rva_vals = h5file["rva_smooth_sat"]
        est_clock_vals = h5file["clk_smooth_sat"]
    else:
        est_rva_vals = h5file["rva_est_sat"]
        est_clock_vals = h5file["clk_est_sat"]

    # Compute the SISE of the position and velocity errors
    errors_rva = est_rva_vals[:, :6] - true_rva_vals[:, :6]

    sise_pos = np.sqrt(
        np.linalg.norm(errors_rva[eval_tidx, :3], axis=1) ** 2
        + (est_clock_vals[eval_tidx, 0] - true_clock_vals[eval_tidx, 0]) ** 2
    )
    sise_vel = (
        np.sqrt(
            np.linalg.norm(errors_rva[eval_tidx, 3:6], axis=1) ** 2
            + (est_clock_vals[eval_tidx, 1] - true_clock_vals[eval_tidx, 1]) ** 2
        )
        * 1000
    )

    rms_pos = np.sqrt(np.mean(sise_pos**2))
    rms_vel = np.sqrt(np.mean(sise_vel**2))

    return sise_pos, sise_vel, rms_pos, rms_vel


def compute_od_errors(
    h5files_filters,
    h5files_smoothers,
    start_ratio,
    end_ratio,
    smooth_selection="best",
    total="mean",
    print_all_cases=True,
):
    # load log files
    n_mc = len(h5files_filters)
    n_iter_save = 10
    print("n_mc:", n_mc)

    n_metrics = 3  # RMS, 95%, 99%

    filter_sise_pos_mat = np.zeros((n_mc + 1, n_iter_save, n_metrics))
    filter_sise_vel_mat = np.zeros((n_mc + 1, n_iter_save, n_metrics))
    filter_sise_pos_cat = np.zeros(0)
    filter_sise_vel_cat = np.zeros(0)
    filter_sise_pos_store = np.zeros((n_mc, n_metrics))
    filter_sise_vel_store = np.zeros((n_mc, n_metrics))
    dz_filter = np.zeros((n_mc + 1, n_iter_save))

    smoother_sise_pos_mat = np.zeros((n_mc + 1, n_iter_save, n_metrics))
    smoother_sise_vel_mat = np.zeros((n_mc + 1, n_iter_save, n_metrics))
    smoother_sise_pos_cat = np.zeros(0)
    smoother_sise_vel_cat = np.zeros(0)
    smoother_sise_pos_store = np.zeros((n_mc, n_metrics))
    smoother_sise_vel_store = np.zeros((n_mc, n_metrics))

    dz_smoother = np.zeros((n_mc + 2, n_iter_save))

    n_mc_valid = 0
    n_iters = [len(h5files_filters[i]) for i in range(n_mc)]
    n_iters.append(1)  # for total

    min_dz_iter = np.zeros(n_mc + 1)

    # Compute Metrics
    for i in range(n_mc):
        n_iter = len(h5files_filters[i])
        exist_smoother = len(h5files_smoothers) > 0
        # if exist_smoother:
        #     # print(h5files_smoothers.keys())
        #     # print(h5files_smoothers[i].keys())
        added_smoother = False

        for j in range(n_iter):
            # Extract Filter Data
            h5file_filter = h5files_filters[i][j]
            tspan = h5file_filter["ts_filter"][:, 0]
            eval_tidx_start = int(len(tspan) * start_ratio)
            eval_tidx_end = int(len(tspan) * end_ratio)
            len_eval = eval_tidx_end - eval_tidx_start

            # if dz_sum_iter does not exist, skip
            if "dz_sum_iter" not in h5file_filter.keys():
                dz_filter[i, j] = 0.0
            else:
                dz_filter[i, j] = h5file_filter["dz_sum_iter"][j, 0]
                if j >= 1:
                    dz_filter[i, j] -= h5file_filter["dz_sum_iter"][0, 0]
            eval_tidx_filter = slice(eval_tidx_start, eval_tidx_end)
            sise_pos, sise_vel, rms_pos, rms_vel = compute_sise_pos_vel(
                h5file_filter, eval_tidx_filter, is_smoother=False
            )

            if rms_pos > 1000:
                # Diverged case - skip
                print(
                    f"MC Run {i+1}, Iteration {j+1} diverged with RMS Position Error {rms_pos:.2f} m. Skipping..."
                )
                continue

            if exist_smoother:
                h5file_smoother = h5files_smoothers[i][j]
                dz_smoother[i, j] = h5file_smoother["dz_sum_iter"][j, 0]
                if j >= 1:
                    dz_smoother[i, j] -= h5file_smoother["dz_sum_iter"][0, 0]
                tspan_smoother = h5file_smoother["ts_filter"][:, 0]
                len_smoother = len(tspan_smoother)
                eval_tidx_smoother = slice(len_smoother - len_eval, len_smoother)
                sise_pos_sm, sise_vel_sm, rms_pos_sm, rms_vel_sm = compute_sise_pos_vel(
                    h5file_smoother, eval_tidx_smoother, is_smoother=True
                )

            # Filter Statistics
            filter_sise_pos_mat[i, j, 0] = rms_pos
            filter_sise_pos_mat[i, j, 1] = np.percentile(sise_pos, 95)
            filter_sise_pos_mat[i, j, 2] = np.percentile(sise_pos, 99.7)
            filter_sise_vel_mat[i, j, 0] = rms_vel
            filter_sise_vel_mat[i, j, 1] = np.percentile(sise_vel, 95)
            filter_sise_vel_mat[i, j, 2] = np.percentile(sise_vel, 99.7)
            # Smoother Statistics
            if exist_smoother:
                smoother_sise_pos_mat[i, j, 0] = rms_pos_sm
                smoother_sise_pos_mat[i, j, 1] = np.percentile(sise_pos_sm, 95)
                smoother_sise_pos_mat[i, j, 2] = np.percentile(sise_pos_sm, 99.7)
                smoother_sise_vel_mat[i, j, 0] = rms_vel_sm
                smoother_sise_vel_mat[i, j, 1] = np.percentile(sise_vel_sm, 95)
                smoother_sise_vel_mat[i, j, 2] = np.percentile(sise_vel_sm, 99.7)

            # register smoother
            register_smoother = False
            if (j == 0) and exist_smoother:
                register_smoother = True
                min_dz_iter[i] = 0
            if j > 0:
                if smooth_selection == "best_dz":
                    # register if dz_filter[i, j] > 0 and dz_smoother[i, j] is minimum
                    if (
                        (dz_filter[i, j] < 0)
                        and (dz_smoother[i, j] < 0)
                        and (dz_smoother[i, j] == np.min(dz_smoother[i, : j + 1]))
                    ):
                        register_smoother = True
                        min_dz_iter[i] = j
                elif smooth_selection == "best":
                    # add if rms_pos_sm is minimum
                    if rms_pos_sm < smoother_sise_pos_mat[i, int(min_dz_iter[i]), 0]:
                        register_smoother = True
                        min_dz_iter[i] = j
                elif smooth_selection == "last":
                    if j == n_iter - 1:
                        register_smoother = True
                        min_dz_iter[i] = j
                elif smooth_selection == "first":
                    # already registered at j==0
                    pass
                else:
                    # directly specify iteration
                    if j == (int(smooth_selection) - 1):
                        register_smoother = True
                        min_dz_iter[i] = j

            if j == 0:
                # For filter, save the first iteration
                filter_sise_pos_cat = np.concatenate((filter_sise_pos_cat, sise_pos))
                filter_sise_vel_cat = np.concatenate((filter_sise_vel_cat, sise_vel))
                filter_sise_pos_store[i, :] = [
                    rms_pos,
                    np.percentile(sise_pos, 95),
                    np.percentile(sise_pos, 99.7),
                ]
                filter_sise_vel_store[i, :] = [
                    rms_vel,
                    np.percentile(sise_vel, 95),
                    np.percentile(sise_vel, 99.7),
                ]
            if register_smoother:
                # For smoother, save the last iteration
                if added_smoother:
                    # if already added previously, remove previous and add new
                    smoother_sise_pos_cat = smoother_sise_pos_cat[:-len_eval]
                    smoother_sise_vel_cat = smoother_sise_vel_cat[:-len_eval]
                smoother_sise_pos_cat = np.concatenate((smoother_sise_pos_cat, sise_pos_sm))
                smoother_sise_vel_cat = np.concatenate((smoother_sise_vel_cat, sise_vel_sm))
                smoother_sise_pos_store[i, :] = [
                    rms_pos_sm,
                    np.percentile(sise_pos_sm, 95),
                    np.percentile(sise_pos_sm, 99.7),
                ]
                smoother_sise_vel_store[i, :] = [
                    rms_vel_sm,
                    np.percentile(sise_vel_sm, 95),
                    np.percentile(sise_vel_sm, 99.7),
                ]
                added_smoother = True

    # Total Statistics
    if total == "cat":
        print("Computing total statistics over all MC runs...")
        filter_sise_pos_mat[n_mc, 0, 0] = np.sqrt(np.mean(filter_sise_pos_cat**2))
        filter_sise_pos_mat[n_mc, 0, 1] = np.percentile(filter_sise_pos_cat, 95)
        filter_sise_pos_mat[n_mc, 0, 2] = np.percentile(filter_sise_pos_cat, 99.7)
        filter_sise_vel_mat[n_mc, 0, 0] = np.sqrt(np.mean(filter_sise_vel_cat**2))
        filter_sise_vel_mat[n_mc, 0, 1] = np.percentile(filter_sise_vel_cat, 95)
        filter_sise_vel_mat[n_mc, 0, 2] = np.percentile(filter_sise_vel_cat, 99.7)
        if exist_smoother:
            print("Computing total statistics for smoother over all MC runs...")
            smoother_sise_pos_mat[n_mc, 0, 0] = np.sqrt(np.mean(smoother_sise_pos_cat**2))
            smoother_sise_pos_mat[n_mc, 0, 1] = np.percentile(smoother_sise_pos_cat, 95)
            smoother_sise_pos_mat[n_mc, 0, 2] = np.percentile(smoother_sise_pos_cat, 99.7)
            smoother_sise_vel_mat[n_mc, 0, 0] = np.sqrt(np.mean(smoother_sise_vel_cat**2))
            smoother_sise_vel_mat[n_mc, 0, 1] = np.percentile(smoother_sise_vel_cat, 95)
            smoother_sise_vel_mat[n_mc, 0, 2] = np.percentile(smoother_sise_vel_cat, 99.7)
    elif total == "mean":
        print("Computing total statistics over mean of MC runs...")
        filter_sise_pos_mat[n_mc, 0, :] = np.mean(filter_sise_pos_store, axis=0)
        filter_sise_vel_mat[n_mc, 0, :] = np.mean(filter_sise_vel_store, axis=0)
        if exist_smoother:
            smoother_sise_pos_mat[n_mc, 0, :] = np.mean(smoother_sise_pos_store, axis=0)
            smoother_sise_vel_mat[n_mc, 0, :] = np.mean(smoother_sise_vel_store, axis=0)
    elif total == "rms":
        print("Computing total statistics over RMS of MC runs...")
        filter_sise_pos_mat[n_mc, 0, :] = np.sqrt(np.mean(filter_sise_pos_store**2, axis=0))
        filter_sise_vel_mat[n_mc, 0, :] = np.sqrt(np.mean(filter_sise_vel_store**2, axis=0))
        if exist_smoother:
            smoother_sise_pos_mat[n_mc, 0, :] = np.sqrt(np.mean(smoother_sise_pos_store**2, axis=0))
            smoother_sise_vel_mat[n_mc, 0, :] = np.sqrt(np.mean(smoother_sise_vel_store**2, axis=0))
    else:
        raise ValueError("Invalid total option. Choose 'cat', 'mean', or 'rms'.")

    # Print statistics
    print("----------------------------------------------------------------------------")
    print("MC    | FS |    Position SISE [m]    |   Velocity SISE [mm/s]   | Residual  ")
    print("Num   |    |  RMS      95%     99.7% |  RMS     95%      99.7%  | Sum       ")
    print("-----------------------------------------------------------------------------")
    for i in range(n_mc + 1):
        mc_str = str(i)
        if not print_all_cases and (i < n_mc):
            continue
        if i == n_mc:
            mc_str = "Total"
        # Print Filter results
        for j in range(n_iters[i]):
            # add star for minimum dz_smoother
            min_dz_symbol = ""
            if exist_smoother and (j == min_dz_iter[i]):
                min_dz_symbol = "*"

            if j == 0:
                print(
                    f"{mc_str:5} | F{j+1:01d} | {filter_sise_pos_mat[i, j, 0]:7.3f} {filter_sise_pos_mat[i, j, 1]:7.3f} {filter_sise_pos_mat[i, j, 2]:7.3f} | {filter_sise_vel_mat[i, j, 0]:7.3f} {filter_sise_vel_mat[i, j, 1]:7.3f} {filter_sise_vel_mat[i, j, 2]:7.3f}  | {dz_filter[i, j]:10.3f}"
                )
            else:
                print(
                    f"      | F{j+1:01d} | {filter_sise_pos_mat[i, j, 0]:7.3f} {filter_sise_pos_mat[i, j, 1]:7.3f} {filter_sise_pos_mat[i, j, 2]:7.3f} | {filter_sise_vel_mat[i, j, 0]:7.3f} {filter_sise_vel_mat[i, j, 1]:7.3f} {filter_sise_vel_mat[i, j, 2]:7.3f}  | {dz_filter[i, j]:10.3f}"
                )
            # Print Smoother results
            if exist_smoother:
                print(
                    f"      | S{j+1:01d} | {smoother_sise_pos_mat[i, j, 0]:7.3f} {smoother_sise_pos_mat[i, j, 1]:7.3f} {smoother_sise_pos_mat[i, j, 2]:7.3f} | {smoother_sise_vel_mat[i, j, 0]:7.3f} {smoother_sise_vel_mat[i, j, 1]:7.3f} {smoother_sise_vel_mat[i, j, 2]:7.3f}  | {dz_smoother[i, j]:10.3f} {min_dz_symbol}"
                )
        print("--------------------------------------------------------------------------------")

    return (
        filter_sise_pos_cat,
        filter_sise_vel_cat,
        smoother_sise_pos_cat,
        smoother_sise_vel_cat,
        filter_sise_pos_mat,
        filter_sise_vel_mat,
        smoother_sise_pos_mat,
        smoother_sise_vel_mat,
    )


#  Plot the error vs covariances
def plot_errors(
    h5file,
    plot_inv=10,
    xlims=None,
    ylims=None,
    use_rtn=True,
    is_smoother=False,
    n_orbits=None,
    fig=None,
    axes=None,
    plot_label="",
    plot_sigma=True,
    sigma_color=None,
):

    #
    tspan = h5file["ts_filter"][::plot_inv, 0]

    true_rva_vals = h5file["rva_true_sat"][::plot_inv, :]
    if is_smoother:
        est_rva_vals = h5file["rva_smooth_sat"][::plot_inv, :]
    else:
        est_rva_vals = h5file["rva_est_sat"][::plot_inv, :]
    rva_sigmas = h5file["rva_sigma_sat"][::plot_inv, :]

    pos_rtn_errs = h5file["/pos_rtn_err"][::plot_inv, :]
    vel_rtn_errs = h5file["/vel_rtn_err"][::plot_inv, :]
    pos_rtn_sigmas = h5file["/pos_rtn_sigma"][::plot_inv, :]
    vel_rtn_sigmas = h5file["/vel_rtn_sigma"][::plot_inv, :]

    true_clock_vals = h5file["clk_true_sat"][::plot_inv, :]
    if is_smoother:
        est_clock_vals = h5file["clk_smooth_sat"][::plot_inv, :]
    else:
        est_clock_vals = h5file["clk_est_sat"][::plot_inv, :]
    clock_sigmas = h5file["clk_sigma_sat"][::plot_inv, :]

    true_srp_vals = h5file["srp_true_sat"][::plot_inv, 0]
    if is_smoother:
        est_srp_vals = h5file["srp_smooth_sat"][::plot_inv, 0]
    else:
        est_srp_vals = h5file["srp_est_sat"][::plot_inv, 0]
    srp_sigmas = h5file["srp_sigma_sat"][::plot_inv, 0]

    # n_meas = h5file["num_meas"]

    if sigma_color is None:
        sigma_color = "blue"

    if ylims is None:
        ylims = {
            "pos": [-1000, 1000],
            "pos_norm": [0, 3000],
            "vel": [-1000, 1000],
            "vel_norm": [0, 3000],
            "clkb": [-1000, 1000],
            "clkd": [-0.1, 0.1],
            "clkdd": [-0.0001, 0.0001],
            "srp": [-1e-2, 1e-2],
        }

    if fig is None or axes is None:
        fig, axes = plt.subplots(3, 4, figsize=(16, 9))

    # Compute the RMS, 95%, 99% of the position, velocity, and clock errors
    sise_pos = np.sqrt(
        np.linalg.norm(est_rva_vals[:, :3] - true_rva_vals[:, :3], axis=1) ** 2
        + (est_clock_vals[:, 0] - true_clock_vals[:, 0]) ** 2
    )
    sise_vel = np.sqrt(
        np.linalg.norm(est_rva_vals[:, 3:6] - true_rva_vals[:, 3:6], axis=1) ** 2
        + (est_clock_vals[:, 1] - true_clock_vals[:, 1]) ** 2
    )
    sise_pos_sigma = np.sqrt(
        np.linalg.norm(rva_sigmas[:, :3], axis=1) ** 2 + (clock_sigmas[:, 0]) ** 2
    )
    sise_vel_sigma = np.sqrt(
        np.linalg.norm(rva_sigmas[:, 3:6], axis=1) ** 2 + (clock_sigmas[:, 1]) ** 2
    )

    if use_rtn:
        posvel_labels = [
            "Radial Position [m]",
            "Tangential Position [m]",
            "Normal Position [m]",
            "SISE Position [m]",
            "Radial Velocity [mm/s]",
            "Tangential Velocity [mm/s]",
            "Normal Velocity [mm/s]",
            "SISE Velocity [mm/s]",
        ]
    else:
        posvel_labels = [
            "Position X [m]",
            "Position Y [m]",
            "Position Z [m]",
            "SISE Position [m]",
            "Velocity X [mm/s]",
            "Velocity Y [mm/s]",
            "Velocity Z [mm/s]",
            "SISE Velocity [mm/s]",
        ]

    ylabels = posvel_labels + [
        "Clock Bias [m]",
        "Clock Drift [m/s]",
        "Clock Drift Rate [m/s^2]",
        "SRP Coeff [m^2/kg]",
    ]

    # Plot the scope of each orbit
    if n_orbits is not None:
        tspan_start = np.zeros(n_orbits)
        for i in range(n_orbits):
            tspan_start[i] = tspan[int(len(tspan) * i / n_orbits)]

    if xlims is None:
        xlims = [tspan[0] / 3600, tspan[-1] / 3600]

    # Fontsizes
    fs = 14
    fs_title = 16
    fs_ticks = 12

    # Position ---------------------------------------------------------------
    for i in range(3):
        if use_rtn:
            errs = pos_rtn_errs[:, i]
            sigmas = pos_rtn_sigmas[:, i]
        else:
            errs = est_rva_vals[:, i] - true_rva_vals[:, i]
            sigmas = rva_sigmas[:, i]
        axes[0, i].plot(tspan / 3600, errs, label=plot_label + " Errors")

        if plot_sigma:
            axes[0, i].fill_between(
                tspan / 3600,
                -3 * sigmas,
                3 * sigmas,
                color=sigma_color,
                alpha=0.1,
                label=plot_label + " 3-sigma",
            )

        # gray background for regions with no measurements
        # no_meas_idx = np.where(n_meas[:, 0] == 0)[0]
        # for idx in no_meas_idx:
        #     axes[0, i].axvspan(tspan[idx] / 3600, tspan[idx] / 3600, color="gray", alpha=0.05, label="No Measurements")

        axes[0, i].set_xlabel("Time [hr]", fontsize=fs)
        axes[0, i].set_ylabel(ylabels[i], fontsize=fs)
        axes[0, i].set_title(ylabels[i], fontsize=fs_title, fontweight="bold")
        axes[0, i].grid(True)
        axes[0, i].legend()
        axes[0, i].set_ylim(ylims["pos"])
        axes[0, i].set_xlim(xlims)
        axes[0, i].tick_params(axis="both", which="major", labelsize=fs_ticks)

    # Norm ---------------------------------------------------------------
    axes[0, 3].plot(
        tspan / 3600,
        sise_pos,
        label=plot_label + " Errors",
    )
    if plot_sigma:
        axes[0, 3].fill_between(
            tspan / 3600,
            0,
            3 * sise_pos_sigma,
            color=sigma_color,
            alpha=0.1,
            label=plot_label + " 3-sigma",
        )
    axes[0, 3].set_xlabel("Time [hr]", fontsize=fs)
    axes[0, 3].set_ylabel(ylabels[3], fontsize=fs)
    axes[0, 3].set_title(ylabels[3], fontsize=fs_title, fontweight="bold")
    axes[0, 3].grid(True)
    axes[0, 3].legend()
    axes[0, 3].set_xlim(xlims)
    axes[0, 3].set_ylim(ylims["pos_norm"])
    axes[0, 3].tick_params(axis="both", which="major", labelsize=fs_ticks)

    if n_orbits is not None:
        for ts in tspan_start:
            # vertical line
            axes[0, 3].axvline(ts / 3600, color="red", linestyle="--")

    # Velocity ---------------------------------------------------------------
    for i in range(3):
        if use_rtn:
            errs = vel_rtn_errs[:, i]
            sigmas = vel_rtn_sigmas[:, i]
        else:
            errs = est_rva_vals[:, i + 3] - true_rva_vals[:, i + 3]
            sigmas = rva_sigmas[:, i + 3]
        axes[1, i].plot(tspan / 3600, errs * 1000, label=plot_label + " Errors")
        if plot_sigma:
            axes[1, i].fill_between(
                tspan / 3600,
                -3 * sigmas * 1000,
                3 * sigmas * 1000,
                color=sigma_color,
                alpha=0.1,
                label=plot_label + " 3-sigma",
            )

        # gray background for regions with no measurements
        # no_meas_idx = np.where(n_meas[:, 0] == 0)[0]
        # for idx in no_meas_idx:
        #     axes[1, i].axvspan(tspan[idx] / 3600, tspan[idx] / 3600, color="gray", alpha=0.05, label="No Measurements")

        axes[1, i].set_xlabel("Time [hr]", fontsize=fs)
        axes[1, i].set_ylabel(ylabels[i + 4], fontsize=fs)
        axes[1, i].set_title(ylabels[i + 4], fontsize=fs_title, fontweight="bold")
        axes[1, i].grid(True)
        axes[1, i].legend()
        axes[1, i].set_ylim(ylims["vel"])
        axes[1, i].set_xlim(xlims)
        axes[1, i].tick_params(axis="both", which="major", labelsize=fs_ticks)

    # Norm ---------------------------------------------------------------
    axes[1, 3].plot(
        tspan / 3600,
        sise_vel * 1000,
        label=plot_label + " Errors",
    )
    if plot_sigma:
        axes[1, 3].fill_between(
            tspan / 3600,
            0,
            3 * sise_vel_sigma * 1000,
            color=sigma_color,
            alpha=0.1,
            label=plot_label + " 3-sigma",
        )
    axes[1, 3].set_xlabel("Time [hr]", fontsize=fs)
    axes[1, 3].set_ylabel(ylabels[7], fontsize=fs)
    axes[1, 3].set_title(ylabels[7], fontsize=fs_title, fontweight="bold")
    axes[1, 3].grid(True)
    axes[1, 3].legend()
    axes[1, 3].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
    axes[1, 3].set_ylim(ylims["vel_norm"])
    axes[1, 3].set_xlim(xlims)
    axes[1, 3].tick_params(axis="both", which="major", labelsize=fs_ticks)

    if n_orbits is not None:
        for ts in tspan_start:
            # vertical line
            axes[1, 3].axvline(ts / 3600, color="red", linestyle="--")

    # Clock ---------------------------------------------------------------
    for i in range(3):
        axes[2, i].plot(
            tspan / 3600, est_clock_vals[:, i] - true_clock_vals[:, i], label=plot_label + " Errors"
        )
        if plot_sigma:
            axes[2, i].fill_between(
                tspan / 3600,
                -3 * clock_sigmas[:, i],
                3 * clock_sigmas[:, i],
                color=sigma_color,
                alpha=0.1,
                label=plot_label + " 3-sigma",
            )

        # gray background for regions with no measurements
        # no_meas_idx = np.where(n_meas[:, 0] == 0)[0]
        # for idx in no_meas_idx:
        #     axes[2, i].axvspan(tspan[idx] / 3600, tspan[idx] / 3600, color="gray", alpha=0.05, label="No Measurements")

        axes[2, i].set_xlabel("Time [hr]", fontsize=fs)
        axes[2, i].set_ylabel(ylabels[i + 8], fontsize=fs)
        axes[2, i].set_title(ylabels[i + 8], fontsize=fs_title, fontweight="bold")
        axes[2, i].grid(True)
        axes[2, i].legend()
        axes[2, i].set_xlim(xlims)
        axes[2, i].tick_params(axis="both", which="major", labelsize=fs_ticks)

        if i == 0:
            axes[2, i].set_ylim(ylims["clkb"])
        elif i == 1:
            axes[2, i].set_ylim(ylims["clkd"])
        else:
            axes[2, i].set_ylim(ylims["clkdd"])

    # SRP ---------------------------------------------------------------
    axes[2, 3].plot(tspan / 3600, est_srp_vals - true_srp_vals, label=plot_label + " Errors")
    if plot_sigma:
        axes[2, 3].fill_between(
            tspan / 3600,
            -3 * srp_sigmas,
            3 * srp_sigmas,
            color=sigma_color,
            alpha=0.1,
            label=plot_label + " 3-sigma",
        )
    axes[2, 3].set_xlabel("Time [hr]", fontsize=fs)
    axes[2, 3].set_ylabel(ylabels[11], fontsize=fs)
    axes[2, 3].set_title(ylabels[11], fontsize=fs_title, fontweight="bold")
    axes[2, 3].grid(True)
    axes[2, 3].legend()
    axes[2, 3].set_xlim(xlims)
    axes[2, 3].set_ylim(ylims["srp"])
    axes[2, 3].set_xlim(xlims)
    axes[2, 3].tick_params(axis="both", which="major", labelsize=fs_ticks)

    # summary
    plt.tight_layout()

    return fig, axes


def plot_sise_errors(
    h5file,
    plot_inv=10,
    xlims=None,
    ylims=None,
    use_rtn=True,
    is_smoother=False,
    n_orbits=None,
    fig=None,
    axes=None,
    plot_label="",
    plot_sigma=True,
    sigma_color=None,
):

    if fig is None or axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3))

    # load vals ----------------------------------------------------------
    tspan = h5file["ts_filter"][::plot_inv, 0]
    true_rva_vals = h5file["rva_true_sat"][::plot_inv, :]
    if is_smoother:
        est_rva_vals = h5file["rva_smooth_sat"][::plot_inv, :]
    else:
        est_rva_vals = h5file["rva_est_sat"][::plot_inv, :]
    rva_sigmas = h5file["rva_sigma_sat"][::plot_inv, :]
    true_clock_vals = h5file["clk_true_sat"][::plot_inv, :]
    if is_smoother:
        est_clock_vals = h5file["clk_smooth_sat"][::plot_inv, :]
    else:
        est_clock_vals = h5file["clk_est_sat"][::plot_inv, :]
    clock_sigmas = h5file["clk_sigma_sat"][::plot_inv, :]

    # Compute the RMS, 95%, 99% of the position, velocity, and clock errors ----
    sise_pos = np.sqrt(
        np.linalg.norm(est_rva_vals[:, :3] - true_rva_vals[:, :3], axis=1) ** 2
        + (est_clock_vals[:, 0] - true_clock_vals[:, 0]) ** 2
    )
    sise_vel = np.sqrt(
        np.linalg.norm(est_rva_vals[:, 3:6] - true_rva_vals[:, 3:6], axis=1) ** 2
        + (est_clock_vals[:, 1] - true_clock_vals[:, 1]) ** 2
    )
    sise_pos_sigma = np.sqrt(
        np.linalg.norm(rva_sigmas[:, :3], axis=1) ** 2 + (clock_sigmas[:, 0]) ** 2
    )
    sise_vel_sigma = np.sqrt(
        np.linalg.norm(rva_sigmas[:, 3:6], axis=1) ** 2 + (clock_sigmas[:, 1]) ** 2
    )

    # Plot the scope of each orbit --------------------------------
    if n_orbits is not None:
        tspan_start = np.zeros(n_orbits)
        for i in range(n_orbits):
            tspan_start[i] = tspan[int(len(tspan) * i / n_orbits)]

    if xlims is None:
        xlims = [tspan[0] / 3600, tspan[-1] / 3600]

    fontsize = 12
    fontsize_ticks = 10

    # Plot SISE position ---------------------------------------------------------------
    axes[0].plot(tspan / 3600, sise_pos, label=plot_label + " Errors")
    if plot_sigma:
        axes[0].fill_between(
            tspan / 3600,
            0,
            3 * sise_pos_sigma,
            color=sigma_color,
            alpha=0.1,
            label=plot_label + " (3-sigma)",
        )
    axes[0].set_xlabel("Time [hr]", fontsize=fontsize)
    axes[0].set_ylabel("SISE Position [m]", fontsize=fontsize)
    axes[0].set_title("SISE Position [m]", fontsize=fontsize)
    axes[0].grid(True)
    axes[0].legend()
    axes[0].set_xlim(xlims)
    axes[0].set_ylim(ylims["pos_norm"])
    # set tick fontsize
    axes[0].tick_params(axis="both", which="major", labelsize=fontsize_ticks)

    if n_orbits is not None:
        for ts in tspan_start:
            # vertical line
            axes[0].axvline(ts / 3600, color="red", linestyle="--")

    # Plot SISE velocity ---------------------------------------------------------------
    axes[1].plot(
        tspan / 3600,
        sise_vel * 1000,
        label=plot_label + " Errors",
    )
    if plot_sigma:
        axes[1].fill_between(
            tspan / 3600,
            0,
            3 * sise_vel_sigma * 1000,
            color=sigma_color,
            alpha=0.1,
            label=plot_label + " (3-sigma)",
        )
    axes[1].set_xlabel("Time [hr]", fontsize=fontsize)
    axes[1].set_ylabel("SISE Velocity [mm/s]", fontsize=fontsize)
    axes[1].set_title("SISE Velocity [mm/s]", fontsize=fontsize)
    axes[1].grid(True)
    axes[1].legend()
    axes[1].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
    axes[1].set_ylim(ylims["vel_norm"])
    axes[1].set_xlim(xlims)
    # set tick fontsize
    axes[1].tick_params(axis="both", which="major", labelsize=fontsize_ticks)

    if n_orbits is not None:
        for ts in tspan_start:
            # vertical line
            axes[1].axvline(ts / 3600, color="red", linestyle="--")

    return fig, axes
