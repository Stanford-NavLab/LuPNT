import numpy as np
import matplotlib.pyplot as plt
import h5py


def compute_od_errors_old(h5files, start_ratio, end_ratio):
    # load logs

    n_mc = len(h5files)
    n_metrics = 3  # RMS, 95%, 99%

    sise_pos_mat = np.zeros((n_mc + 1, n_metrics))
    sise_vel_mat = np.zeros((n_mc + 1, n_metrics))
    sise_pos_cat = np.zeros(0)
    sise_vel_cat = np.zeros(0)

    n_mc_valid = 0

    # Compute Metrics
    for i, h5file in enumerate(h5files):
        tspan = h5file["ts_filter"][:, 0]
        true_rva_vals = h5file["rva_true_sat"]
        est_rva_vals = h5file["rva_est_sat"]

        true_clock_vals = h5file["clk_true_sat"]
        est_clock_vals = h5file["clk_est_sat"]

        true_srp_vals = h5file["srp_true_sat"][:, 0]
        est_srp_vals = h5file["srp_est_sat"][:, 0]

        # Compute the RMS, 95%, 99% of the position, velocity, and clock errors
        errors_rva = est_rva_vals[:, :6] - true_rva_vals[:, :6]
        eval_tidx_start = int(len(tspan) * start_ratio)
        eval_tidx_end = int(len(tspan) * end_ratio)
        eval_tidx = slice(eval_tidx_start, eval_tidx_end)

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
        if rms_pos > 1000 or rms_vel > 1000:
            print(f"Warning: Large SISE detected in MC {i}: Pos {rms_pos} m, Vel {rms_vel} mm/s")
            continue
        if rms_pos < 1e-2 or rms_vel < 1e-2:
            print(f"Warning: Small SISE detected in MC {i}: Pos {rms_pos} m, Vel {rms_vel} mm/s")
            continue

        sise_pos_cat = np.concatenate((sise_pos_cat, sise_pos))
        sise_vel_cat = np.concatenate((sise_vel_cat, sise_vel))

        sise_pos_mat[i, 0] = rms_pos
        sise_pos_mat[i, 1] = np.percentile(sise_pos, 95)
        sise_pos_mat[i, 2] = np.percentile(sise_pos, 99)
        sise_vel_mat[i, 0] = rms_vel
        sise_vel_mat[i, 1] = np.percentile(sise_vel, 95)
        sise_vel_mat[i, 2] = np.percentile(sise_vel, 99)

        n_mc_valid += 1

    # Total Statistics
    sise_pos_mat[n_mc_valid, 0] = np.sqrt(np.mean(sise_pos_cat**2))
    sise_pos_mat[n_mc_valid, 1] = np.percentile(sise_pos_cat, 95)
    sise_pos_mat[n_mc_valid, 2] = np.percentile(sise_pos_cat, 99)
    sise_vel_mat[n_mc_valid, 0] = np.sqrt(np.mean(sise_vel_cat**2))
    sise_vel_mat[n_mc_valid, 1] = np.percentile(sise_vel_cat, 95)
    sise_vel_mat[n_mc_valid, 2] = np.percentile(sise_vel_cat, 99)

    # Print statistics
    print("-----------------------------------------------------------")
    print("MC    |     Position SISE [m]   |   Velocity SISE [mm/s]     ")
    print("Num   | RMS      95%     99%  |  RMS     95%      99%      ")
    print("-----------------------------------------------------------")
    for i in range(n_mc_valid + 1):
        mc_str = str(i)
        if i == n_mc_valid:
            mc_str = "Total"
            print("-----------------------------------------------------------")
        print(
            f"{mc_str:5} | {sise_pos_mat[i, 0]:7.3f} {sise_pos_mat[i, 1]:7.3f} {sise_pos_mat[i, 2]:7.3f} | {sise_vel_mat[i, 0]:7.3f} {sise_vel_mat[i, 1]:7.3f} {sise_vel_mat[i, 2]:7.3f}"
        )
    print("-----------------------------------------------------------")

    return sise_pos_mat, sise_vel_mat


#  Plot the error vs covariances
def plot_errors_old(h5file, plot_inv=10, ylims=None, use_rtn=True, n_orbits=6):

    #
    tspan = h5file["ts_filter"][::plot_inv, 0]
    true_rva_vals = h5file["rva_true_sat"][::plot_inv, :]
    est_rva_vals = h5file["rva_est_sat"][::plot_inv, :]
    rva_sigmas = h5file["rva_sigma_sat"][::plot_inv, :]

    pos_rtn_errs = h5file["/pos_rtn_err"][::plot_inv, :]
    vel_rtn_errs = h5file["/vel_rtn_err"][::plot_inv, :]
    pos_rtn_sigmas = h5file["/pos_rtn_sigma"][::plot_inv, :]
    vel_rtn_sigmas = h5file["/vel_rtn_sigma"][::plot_inv, :]

    true_clock_vals = h5file["clk_true_sat"][::plot_inv, :]
    est_clock_vals = h5file["clk_est_sat"][::plot_inv, :]
    clock_sigmas = h5file["clk_sigma_sat"][::plot_inv, :]

    true_srp_vals = h5file["srp_true_sat"][::plot_inv, 0]
    est_srp_vals = h5file["srp_est_sat"][::plot_inv, 0]
    srp_sigmas = h5file["srp_sigma_sat"][::plot_inv, 0]

    # n_meas = h5file["num_meas"]

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
            "Position R [m]",
            "Position T [m]",
            "Position N [m]",
            "SISE Position [m]",
            "Velocity R [mm/s]",
            "Velocity T [mm/s]",
            "Velocity N [mm/s]",
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
        "SRP Coeffs [m/s^2]",
    ]

    # Plot the scope of each orbit
    tspan_start = np.zeros(n_orbits)
    for i in range(n_orbits):
        tspan_start[i] = tspan[int(len(tspan) * i / n_orbits)]

    # Position ---------------------------------------------------------------
    for i in range(3):
        if use_rtn:
            errs = pos_rtn_errs[:, i]
            sigmas = pos_rtn_sigmas[:, i]
        else:
            errs = est_rva_vals[:, i] - true_rva_vals[:, i]
            sigmas = rva_sigmas[:, i]
        axes[0, i].plot(tspan / 3600, errs, label="Errors")
        axes[0, i].fill_between(
            tspan / 3600,
            -3 * sigmas,
            3 * sigmas,
            color="blue",
            alpha=0.1,
            label="3-sigma",
        )

        # gray background for regions with no measurements
        # no_meas_idx = np.where(n_meas[:, 0] == 0)[0]
        # for idx in no_meas_idx:
        #     axes[0, i].axvspan(tspan[idx] / 3600, tspan[idx] / 3600, color="gray", alpha=0.05, label="No Measurements")

        axes[0, i].set_xlabel("Time [hr]")
        axes[0, i].set_ylabel(ylabels[i])
        axes[0, i].set_title(ylabels[i])
        axes[0, i].grid(True)
        axes[0, i].legend()
        axes[0, i].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
        axes[0, i].set_ylim(ylims["pos"])

    # Norm ---------------------------------------------------------------
    axes[0, 3].plot(
        tspan / 3600,
        sise_pos,
        label="Errors",
    )
    axes[0, 3].fill_between(
        tspan / 3600,
        0,
        3 * sise_pos_sigma,
        color="blue",
        alpha=0.1,
        label="3-sigma",
    )
    axes[0, 3].set_xlabel("Time [hr]")
    axes[0, 3].set_ylabel(ylabels[3])
    axes[0, 3].set_title(ylabels[3])
    axes[0, 3].grid(True)
    axes[0, 3].legend()
    axes[0, 3].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
    axes[0, 3].set_ylim(ylims["pos_norm"])
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
        axes[1, i].plot(tspan / 3600, errs * 1000, label="Errors")
        axes[1, i].fill_between(
            tspan / 3600,
            -3 * sigmas * 1000,
            3 * sigmas * 1000,
            color="blue",
            alpha=0.1,
            label="3-sigma",
        )

        # gray background for regions with no measurements
        # no_meas_idx = np.where(n_meas[:, 0] == 0)[0]
        # for idx in no_meas_idx:
        #     axes[1, i].axvspan(tspan[idx] / 3600, tspan[idx] / 3600, color="gray", alpha=0.05, label="No Measurements")

        axes[1, i].set_xlabel("Time [hr]")
        axes[1, i].set_ylabel(ylabels[i + 4])
        axes[1, i].set_title(ylabels[i + 4])
        axes[1, i].grid(True)
        axes[1, i].legend()
        axes[1, i].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
        axes[1, i].set_ylim(ylims["vel"])

    # Norm ---------------------------------------------------------------
    axes[1, 3].plot(
        tspan / 3600,
        sise_vel * 1000,
        label="Errors",
    )
    axes[1, 3].fill_between(
        tspan / 3600,
        0,
        3 * sise_vel_sigma * 1000,
        color="blue",
        alpha=0.1,
        label="3-sigma",
    )
    axes[1, 3].set_xlabel("Time [hr]")
    axes[1, 3].set_ylabel(ylabels[7])
    axes[1, 3].set_title(ylabels[7])
    axes[1, 3].grid(True)
    axes[1, 3].legend()
    axes[1, 3].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
    axes[1, 3].set_ylim(ylims["vel_norm"])
    for ts in tspan_start:
        # vertical line
        axes[1, 3].axvline(ts / 3600, color="red", linestyle="--")

    # Clock ---------------------------------------------------------------
    for i in range(3):
        axes[2, i].plot(tspan / 3600, est_clock_vals[:, i] - true_clock_vals[:, i], label="Errors")
        axes[2, i].fill_between(
            tspan / 3600,
            -3 * clock_sigmas[:, i],
            3 * clock_sigmas[:, i],
            color="blue",
            alpha=0.1,
            label="3-sigma",
        )

        # gray background for regions with no measurements
        # no_meas_idx = np.where(n_meas[:, 0] == 0)[0]
        # for idx in no_meas_idx:
        #     axes[2, i].axvspan(tspan[idx] / 3600, tspan[idx] / 3600, color="gray", alpha=0.05, label="No Measurements")

        axes[2, i].set_xlabel("Time [hr]")
        axes[2, i].set_ylabel(ylabels[i + 8])
        axes[2, i].set_title(ylabels[i + 8])
        axes[2, i].grid(True)
        axes[2, i].legend()
        axes[2, i].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
        if i == 0:
            axes[2, i].set_ylim(ylims["clkb"])
        elif i == 1:
            axes[2, i].set_ylim(ylims["clkd"])
        else:
            axes[2, i].set_ylim(ylims["clkdd"])

    # SRP ---------------------------------------------------------------
    axes[2, 3].plot(tspan / 3600, est_srp_vals - true_srp_vals, label="Errors")
    axes[2, 3].fill_between(
        tspan / 3600,
        -3 * srp_sigmas,
        3 * srp_sigmas,
        color="blue",
        alpha=0.1,
        label="3-sigma",
    )
    axes[2, 3].set_xlabel("Time [hr]")
    axes[2, 3].set_ylabel(ylabels[11])
    axes[2, 3].set_title(ylabels[11])
    axes[2, 3].grid(True)
    axes[2, 3].legend()
    axes[2, 3].set_xlim([tspan[0] / 3600, tspan[-1] / 3600])
    axes[2, 3].set_ylim(ylims["srp"])

    # summary
    plt.tight_layout()
    plt.show()
