import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from pylupnt.numerics import metrics


def plot_ephemeris_error(dfs, diff_norm=False, inv=1, filename=None):

    # plot for different constelaltions
    colors = ["blue", "orange", "green", "cyan", "magenta", "yellow"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fs = 16  # fontsize
    fs_legend = 14
    fs_ticks = 14
    # set tick param size
    for ax in axes:
        ax.tick_params(axis="both", labelsize=fs_ticks)  # 'x', 'y', or 'both'  # font size

    gnss_consts = list(dfs.keys())
    for i, gnss_const in enumerate(gnss_consts):
        df = dfs[gnss_const][1]  # get dataframe
        pos_tx = df["pos_tx"].values[::inv]
        vel_tx = df["vel_tx"].values[::inv]
        pos_rx = df["pos_rx"].values[::inv]
        clock_tx = df["clockbias_tx"].values[::inv]
        pos_tx_ephem = df["pos_tx_ephem"].values[::inv]
        vel_tx_ephem = df["vel_tx_ephem"].values[::inv]
        clock_tx_ephem = df["clockbias_ephem"].values[::inv]
        pco_x = df["pco_ecef_x_m"].values[::inv]
        pco_y = df["pco_ecef_y_m"].values[::inv]
        pco_z = df["pco_ecef_z_m"].values[::inv]
        pco_ecef_m = np.vstack((pco_x, pco_y, pco_z)).T  # shape (N, 3)

        clock_tx_diff_median = df["clock_bias_tx_median"].values[
            ::inv
        ]  # median clock bias in seconds

        # plot differences
        if diff_norm:
            # simply use norm
            pos_tx_diff = (
                np.array(
                    [
                        np.linalg.norm(pos_tx[i] + pco_ecef_m[i] / 1000)
                        - np.linalg.norm(pos_tx_ephem[i])
                        for i in range(len(pos_tx))
                    ]
                )
                * 1000
            )
        else:
            pos_tx_diff = (
                np.array(
                    [
                        np.linalg.norm(pos_tx[i] + pco_ecef_m[i] / 1000 - pos_rx[i])
                        - np.linalg.norm(pos_tx_ephem[i] - pos_rx[i])
                        for i in range(len(pos_tx))
                    ]
                )
                * 1000
            )
        C_ms = 299792458.0  # speed of light in m/s
        clock_tx_diff = (clock_tx - clock_tx_ephem + clock_tx_diff_median) * C_ms
        ure_tx_diff = pos_tx_diff + clock_tx_diff

        pos_tx_diff_without_outliers = pos_tx_diff[np.abs(pos_tx_diff) < 10]
        ure_tx_diff_without_outliers = ure_tx_diff[np.abs(ure_tx_diff) < 10]
        clock_tx_diff_without_outliers = clock_tx_diff[np.abs(clock_tx_diff) < 10]

        print(
            f"[Position] GNSS: {gnss_const}  mean: ",
            np.mean(pos_tx_diff_without_outliers),
            " m,  std: ",
            np.std(pos_tx_diff_without_outliers),
            " m",
        )
        print(
            f"[Clock] GNSS: {gnss_const}  mean: ",
            np.mean(clock_tx_diff_without_outliers),
            " m,  std: ",
            np.std(clock_tx_diff_without_outliers),
            " m",
        )
        print(
            f"[Total] GNSS: {gnss_const}  mean: ",
            np.mean(ure_tx_diff_without_outliers),
            " m,  std: ",
            np.std(ure_tx_diff_without_outliers),
            " m",
        )

        bins_pos = np.linspace(-5, 5, 100)  # remove outliers
        axes[0].hist(
            pos_tx_diff, bins=bins_pos, color=colors[i], alpha=0.7, label=gnss_const, density=False
        )
        axes[0].set_title("Transmitter Position Difference", fontsize=fs)
        axes[0].set_xlabel("Position Difference (m)", fontsize=fs)
        axes[0].set_ylabel("Count", fontsize=fs)
        axes[0].grid(True)
        axes[0].legend(fontsize=fs_legend)

        bins_clk = np.linspace(-5, 5, 100)  # remove outliers
        axes[1].hist(
            clock_tx_diff,
            bins=bins_clk,
            color=colors[i],
            alpha=0.7,
            label=gnss_const,
            density=False,
        )
        axes[1].set_title("Transmitter Clock Bias Difference", fontsize=fs)
        axes[1].set_xlabel("Clock Bias Difference (m)", fontsize=fs)
        axes[1].set_ylabel("Count", fontsize=fs)
        axes[1].grid(True)
        axes[1].legend(fontsize=fs_legend)

        bins_vel = np.linspace(-5, 5, 100)  # remove outliers
        axes[2].hist(
            ure_tx_diff, bins=bins_vel, color=colors[i], alpha=0.7, label=gnss_const, density=False
        )
        axes[2].set_title("Transmitter URE Difference", fontsize=fs)
        axes[2].set_xlabel("URE Difference (m)", fontsize=fs)
        axes[2].set_ylabel("Count", fontsize=fs)
        axes[2].grid(True)
        axes[2].legend(fontsize=fs_legend)

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300)
    plt.show()


def generate_timestep_df(df_signals, gnss_const, signal_families, tidx_inv=60):
    # Generate df that includes time step information
    df_timestep = {}

    if tidx_inv == 1:
        # not allowed
        raise ValueError("tidx_inv cannot be 1, must be greater than 1 to perform subsampling.")

    # Initialize df ----------------------------------------------------------------
    for signal_idx in signal_families:
        # construct dateframe
        df = df_signals[gnss_const][signal_idx]
        tidxs_load = df["tidx"].values
        max_tidx = tidxs_load[-1]  # assuming tidx is sorted
        tidxs = np.arange(1, max_tidx + 1, tidx_inv)
        lent = len(tidxs)

        total_meas = len(df_signals[gnss_const][signal_idx])
        df_timestep[signal_idx] = {
            "lent": lent,
            "tidxs": tidxs,
            "num_sats": np.zeros(lent, dtype=int),
            "pr_noise": np.zeros(total_meas, dtype=float),
            "pr_iono_error": np.zeros(total_meas, dtype=float),
            "pr_ure_error": np.zeros(total_meas, dtype=float),
            "tdcp_noise": np.zeros(total_meas, dtype=float),
            "tdcp_ure_error": np.zeros(total_meas, dtype=float),
            "tdcp_iono_error": np.zeros(total_meas, dtype=float),
            "tdcp_velocity": np.zeros(total_meas, dtype=float),
            "tdcp_relpos": np.zeros(total_meas, dtype=float),
            "min_alt_pr": np.zeros(total_meas, dtype=float),
            "min_alt_tdcp": np.zeros(total_meas, dtype=float),
            "cn0": np.zeros(total_meas, dtype=float),
            "tracked_prns": [[] for _ in range(lent)],
            "tdcp_prns": [[] for _ in range(lent)],
            "prn_rows_tidx": [{} for _ in range(lent)],
            "prn_measid_tidx": [{} for _ in range(lent)],
        }

    df_timestep["L1L5"] = {
        "num_sats": np.zeros(lent, dtype=int),
        "delays_iono": np.zeros(total_meas, dtype=float),
        "min_alt": np.zeros(total_meas, dtype=float),
        "pr_noise": np.zeros(total_meas, dtype=float),
        "cn0": np.zeros(total_meas, dtype=float),
        "ure_error": np.zeros(total_meas, dtype=float),
    }

    # Compute Delays for Each Frequency ----------------------------------------------
    for signal_idx in signal_families:
        # timesteps
        df = df_signals[gnss_const][signal_idx]
        tidxs_load = df["tidx"].values
        max_tidx = tidxs_load[-1]  # assuming tidx is sorted
        tidxs = np.arange(1, max_tidx + 1, tidx_inv)
        lent = len(tidxs)

        print("----------------------------------------")
        print("Processing constellation:", gnss_const, "signal index:", signal_idx)
        print("-----------------------------------------")

        # Extract relevant columns with subsampling
        total_delay_m = df["total_delay_m"].values
        sigma_range = df["sigma_range"].values * 1000  # convert to meters
        sigma_carrier = df["sigma_carrier"].values * 1000  # convert to meters
        min_alts = df["min_alt"].values
        pos_tx = df["pos_tx"].values
        vel_tx = df["vel_tx"].values
        pos_rx = df["pos_rx"].values
        vel_rx = df["vel_rx"].values
        clock_tx = df["clockbias_tx"].values
        pos_tx_ephem = df["pos_tx_ephem"].values
        vel_tx_ephem = df["vel_tx_ephem"].values
        pco_x = df["pco_ecef_x_m"].values
        pco_y = df["pco_ecef_y_m"].values
        pco_z = df["pco_ecef_z_m"].values
        pco_ecef = np.vstack((pco_x, pco_y, pco_z)).T  # shape (N, 3)
        clock_tx_ephem = df["clockbias_ephem"].values
        cn0 = df["cn0"].values

        prn = df["prn"].values

        # Compute ephemeris errors
        pos_tx_diff = (
            np.array(
                [
                    np.linalg.norm(pos_tx[i] + pco_ecef[i] / 1000 - pos_rx[i])
                    - np.linalg.norm(pos_tx_ephem[i] - pos_rx[i])
                    for i in range(len(pos_tx))
                ]
            )
            * 1000
        )
        C_ms = 299792458.0  # speed of light in m/s
        clock_tx_diff = (clock_tx - clock_tx_ephem) * C_ms
        ure_tx_diff = pos_tx_diff + clock_tx_diff
        print(
            "ure_tx_diff stats (m): mean =",
            np.median(ure_tx_diff),
            ", 68% =",
            np.percentile(np.abs(ure_tx_diff), 68),
        )

        curr_row = 0
        num_meas = 0
        num_meas_tdcp = 0

        prev_tracked_prns = set()
        prev_prn_rows = {}

        for k, tidx in tqdm(
            enumerate(tidxs), total=lent, desc=f"GNSS: {gnss_const}, Signal: {signal_idx}"
        ):
            # find the rows corresponding to the current tidx
            start_row = curr_row
            while True:
                if curr_row >= len(df):
                    # reached the end of dataframe
                    break

                row_tidx = df["tidx"].values[curr_row]
                if row_tidx == tidx:
                    # timestep match
                    df_timestep[signal_idx]["tracked_prns"][k].append(prn[curr_row])
                    df_timestep[signal_idx]["prn_rows_tidx"][k][prn[curr_row]] = curr_row
                    df_timestep[signal_idx]["pr_ure_error"][num_meas] = ure_tx_diff[curr_row]
                    df_timestep[signal_idx]["pr_noise"][num_meas] = sigma_range[curr_row]
                    df_timestep[signal_idx]["pr_iono_error"][num_meas] = total_delay_m[curr_row]
                    df_timestep[signal_idx]["min_alt_pr"][num_meas] = min_alts[curr_row]
                    df_timestep[signal_idx]["cn0"][num_meas] = cn0[curr_row]

                    df_timestep[signal_idx]["prn_measid_tidx"][k][prn[curr_row]] = num_meas

                    # Check if the prn is also tracked on previous timestep
                    if k > 0:
                        # If TDCP is observable, add delay here
                        if prn[curr_row] in prev_tracked_prns:
                            prev_row = prev_prn_rows.get(prn[curr_row], [])
                            df_timestep[signal_idx]["tdcp_prns"][k].append(prn[curr_row])

                            curr_noise = sigma_carrier[curr_row]
                            prev_noise = sigma_carrier[prev_row]
                            df_timestep[signal_idx]["tdcp_noise"][num_meas_tdcp] = np.sqrt(
                                curr_noise**2 + prev_noise**2
                            )

                            curr_ure_error = ure_tx_diff[curr_row]
                            prev_ure_error = ure_tx_diff[prev_row]
                            df_timestep[signal_idx]["tdcp_ure_error"][num_meas_tdcp] = (
                                curr_ure_error - prev_ure_error
                            )

                            curr_delay_iono = total_delay_m[curr_row]
                            prev_delay_iono = total_delay_m[prev_row]
                            df_timestep[signal_idx]["tdcp_iono_error"][num_meas_tdcp] = (
                                curr_delay_iono - prev_delay_iono
                            )

                            prev_vel_rx = vel_rx[prev_row]
                            df_timestep[signal_idx]["tdcp_velocity"][num_meas_tdcp] = (
                                np.linalg.norm(prev_vel_rx)
                            )

                            curr_pos_tx = pos_tx[curr_row]
                            prev_pos_tx = pos_tx[prev_row]
                            curr_relpos = np.linalg.norm(curr_pos_tx - pos_rx[curr_row])
                            prev_relpos = np.linalg.norm(prev_pos_tx - pos_rx[prev_row])
                            df_timestep[signal_idx]["tdcp_relpos"][num_meas_tdcp] = (
                                curr_relpos - prev_relpos
                            )

                            df_timestep[signal_idx]["min_alt_tdcp"][num_meas_tdcp] = min_alts[
                                curr_row
                            ]

                            num_meas_tdcp += 1

                    curr_row += 1
                    num_meas += 1
                elif row_tidx > tidx:
                    # has entered next timestep
                    break
                else:
                    # less than current timestep, keep searching
                    if row_tidx == tidx - 1:
                        prev_tracked_prns.add(prn[curr_row])
                        prev_prn_rows[prn[curr_row]] = curr_row
                    curr_row += 1

            # Store number of tracked satellites
            df_timestep[signal_idx]["num_sats"][k] = len(df_timestep[signal_idx]["tracked_prns"][k])

            # reset for next timestep
            prev_tracked_prns = set()
            prev_prn_rows = {}

        # Remove redundant entries in tdcp_prns and resize noise/error arrays
        df_timestep[signal_idx]["pr_noise"] = df_timestep[signal_idx]["pr_noise"][:num_meas]
        df_timestep[signal_idx]["pr_iono_error"] = df_timestep[signal_idx]["pr_iono_error"][
            :num_meas
        ]
        df_timestep[signal_idx]["min_alt_pr"] = df_timestep[signal_idx]["min_alt_pr"][:num_meas]
        df_timestep[signal_idx]["cn0"] = df_timestep[signal_idx]["cn0"][:num_meas]
        df_timestep[signal_idx]["pr_ure_error"] = df_timestep[signal_idx]["pr_ure_error"][:num_meas]
        df_timestep[signal_idx]["tdcp_noise"] = df_timestep[signal_idx]["tdcp_noise"][
            :num_meas_tdcp
        ]
        df_timestep[signal_idx]["tdcp_ure_error"] = df_timestep[signal_idx]["tdcp_ure_error"][
            :num_meas_tdcp
        ]
        df_timestep[signal_idx]["tdcp_iono_error"] = df_timestep[signal_idx]["tdcp_iono_error"][
            :num_meas_tdcp
        ]
        df_timestep[signal_idx]["tdcp_velocity"] = df_timestep[signal_idx]["tdcp_velocity"][
            :num_meas_tdcp
        ]
        df_timestep[signal_idx]["tdcp_relpos"] = df_timestep[signal_idx]["tdcp_relpos"][
            :num_meas_tdcp
        ]
        df_timestep[signal_idx]["min_alt_tdcp"] = df_timestep[signal_idx]["min_alt_tdcp"][
            :num_meas_tdcp
        ]

    # Ionofree combinations for L1L5 ------------------------------------------------
    lent_L1 = df_timestep[1]["lent"]
    lent_L5 = df_timestep[5]["lent"]
    lent_common = min(lent_L1, lent_L5)
    df_timestep["L1L5"]["tidxs"] = np.zeros(lent_common, dtype=int)
    df_timestep["L1L5"]["num_sats"] = np.zeros(lent_common, dtype=int)

    f_L1 = 1575.42  # L1 frequency in Hz
    f_L5 = 1176.45  # L5 frequency in Hz
    alpha_L1 = f_L1**2 / (f_L1**2 - f_L5**2)
    alpha_L5 = f_L5**2 / (f_L1**2 - f_L5**2)

    num_meas_ionofree = 0

    for k in tqdm(range(lent_common), desc="Computing L1L5 common sats"):
        prns_L1 = set(df_timestep[1]["tracked_prns"][k])
        prns_L5 = set(df_timestep[5]["tracked_prns"][k])
        prns_common = prns_L1.intersection(prns_L5)
        df_timestep["L1L5"]["num_sats"][k] = len(prns_common)
        if lent_L1 > lent_L5:
            # L5 shorter, copy from L5
            df_timestep["L1L5"]["tidxs"][k] = df_timestep[5]["tidxs"][k]
        else:
            # L1 shorter, copy from L1
            df_timestep["L1L5"]["tidxs"][k] = df_timestep[1]["tidxs"][k]

        # Compute ionofree delays and min altitudes
        for prn in prns_common:
            row_L1 = df_timestep[1]["prn_rows_tidx"][k][prn]
            row_L5 = df_timestep[5]["prn_rows_tidx"][k][prn]

            delay_iono_L1 = df_signals[gnss_const][1]["total_delay_m"][row_L1]
            delay_iono_L5 = df_signals[gnss_const][5]["total_delay_m"][row_L5]
            pr_noise_L1 = df_signals[gnss_const][1]["sigma_range"][row_L1] * 1000  # in meters
            pr_noise_L5 = df_signals[gnss_const][5]["sigma_range"][row_L5] * 1000  # in meters

            # ionofree delay
            ionofree_delay = alpha_L1 * delay_iono_L1 - alpha_L5 * delay_iono_L5
            df_timestep["L1L5"]["delays_iono"][
                num_meas_ionofree
            ] = ionofree_delay  # store using L1 row index

            # ionofree noise
            ionofree_noise = np.sqrt((alpha_L1 * pr_noise_L1) ** 2 + (alpha_L5 * pr_noise_L5) ** 2)
            df_timestep["L1L5"]["pr_noise"][
                num_meas_ionofree
            ] = ionofree_noise  # store using L1 row index

            # min altitude (take from L1)
            min_alt = df_signals[gnss_const][1]["min_alt"][row_L1]
            df_timestep["L1L5"]["min_alt"][num_meas_ionofree] = min_alt  # store using L1 row index

            # cn0 (take from L1)
            cn0 = df_signals[gnss_const][1]["cn0"][row_L1]
            df_timestep["L1L5"]["cn0"][num_meas_ionofree] = cn0  # store using L1 row index

            # ephemeris error
            measid_L1 = df_timestep[1]["prn_measid_tidx"][k][prn]
            measid_L5 = df_timestep[5]["prn_measid_tidx"][k][prn]
            ure_error_L1 = df_timestep[1]["pr_ure_error"][measid_L1]
            ure_error_L5 = df_timestep[5]["pr_ure_error"][measid_L5]
            ionofree_ure_error = alpha_L1 * ure_error_L1 - alpha_L5 * ure_error_L5
            df_timestep["L1L5"]["ure_error"][
                num_meas_ionofree
            ] = ionofree_ure_error  # store using L1 row index

            num_meas_ionofree += 1

    # remove redundant entries in L1L5
    df_timestep["L1L5"]["delays_iono"] = df_timestep["L1L5"]["delays_iono"][:num_meas_ionofree]
    df_timestep["L1L5"]["min_alt"] = df_timestep["L1L5"]["min_alt"][:num_meas_ionofree]
    df_timestep["L1L5"]["pr_noise"] = df_timestep["L1L5"]["pr_noise"][:num_meas_ionofree]
    df_timestep["L1L5"]["cn0"] = df_timestep["L1L5"]["cn0"][:num_meas_ionofree]
    df_timestep["L1L5"]["ure_error"] = df_timestep["L1L5"]["ure_error"][:num_meas_ionofree]

    print(
        "Ionofree URE STATS (m): mean =",
        np.median(df_timestep["L1L5"]["ure_error"]),
        ", 68% =",
        np.percentile(np.abs(df_timestep["L1L5"]["ure_error"]), 68),
    )

    return df_timestep


def plot_tracked_sats(df_timestep, gnss_consts, filename=None):
    # For L1, L5, and L1L5 ionofree
    fig, axes = plt.subplots(2, 3, figsize=(16, 6))

    signal_families = [1, 5, "L1L5"]
    label_signals = ["L1/E1", "L5/E5a", "Ionofree"]
    colors = ["blue", "orange", "green", "black"]

    fs_ticks = 14
    fs_labels = 14

    for j, signal_idx in enumerate(signal_families):
        total_sat_nums = {}

        for i, gnss_const in enumerate(gnss_consts):
            num_sats = df_timestep[gnss_const][signal_idx]["num_sats"]
            tidxs = df_timestep[gnss_const][signal_idx]["tidxs"]

            # store common tidxs
            for k, tidx in enumerate(tidxs):
                # if not in total_sat_nums, initialize
                if tidx not in total_sat_nums:
                    total_sat_nums[tidx] = num_sats[k]
                else:
                    total_sat_nums[tidx] += num_sats[k]

            axes[0, j].plot(tidxs / 3600, num_sats, label=gnss_const, color=colors[i], alpha=0.7)
            axes[0, j].set_xlabel("Time [hr]", fontsize=16)
            axes[0, j].set_ylabel("Counts", fontsize=16)
            axes[0, j].set_title(f"{label_signals[j]}", fontsize=16, fontweight="bold")
            axes[0, j].grid(True)
            axes[0, j].legend(fontsize=fs_labels)
            axes[0, j].tick_params(axis="both", labelsize=fs_ticks)
            axes[0, j].set_ylim([0, 15])
            axes[0, j].set_yticks(np.arange(0, 16, 2))

        # construct total tracked sats -------------------------------------------------
        total_tidxs = sorted(total_sat_nums.keys())
        total_sats = [total_sat_nums[tidx] for tidx in total_tidxs]
        axes[1, j].plot(np.array(total_tidxs) / 3600, total_sats, label="Total", color=colors[-1])
        axes[1, j].set_xlabel("Time [hr]", fontsize=16)
        axes[1, j].set_ylabel("Counts", fontsize=16)
        axes[1, j].set_title(f"{label_signals[j]} (Total)", fontsize=16, fontweight="bold")
        axes[1, j].grid(True)
        axes[1, j].legend(fontsize=fs_labels)
        axes[1, j].tick_params(axis="both", labelsize=fs_ticks)
        axes[1, j].set_ylim([0, 20])
        axes[1, j].set_yticks(np.arange(0, 21, 2))

        print(f"Average number of tracked sats for {label_signals[j]}: ", np.mean(total_sats))
        print(f"Max number of tracked sats for {label_signals[j]}: ", np.max(total_sats))
        print(f"Standard deviation of tracked sats for {label_signals[j]}: ", np.std(total_sats))

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300)
    plt.show()


def plot_tdcp_errors(df_timestep, gnss_consts, signal_family, plot_inv=10, filename=None):
    # Plot TDCP errors for different constellations and signal families
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))
    fs = 16
    fs_legend = 14
    fs_ticks = 14
    fs_title = 16
    xlims = [0, 6]  # time in hours

    ylabels = [
        "TDCP Receiver Noise (mm)",
        "TDCP Ephemeris Error (mm)",
        "TDCP Ionospheric Delay (mm)",
    ]
    metrics_plot = ["tdcp_noise", "tdcp_ure_error", "tdcp_iono_error"]
    colors = ["blue", "orange", "green", "purple"]
    range_metrics = [[0, 10], [-10, 10], [-10, 10]]  # in mm
    alt_ranges = [[0, 1000], [1000, 5000], [5000, 20000]]  # in km

    for i, alt_ranges in enumerate(alt_ranges):
        plot_vals = {}
        for j, gnss_const in enumerate(gnss_consts):
            for k, metric in enumerate(metrics_plot):
                # add to dictionary
                alt_ranges_km = df_timestep[gnss_const][signal_family]["min_alt_tdcp"][::plot_inv]
                alt_mask = (alt_ranges_km >= alt_ranges[0]) & (alt_ranges_km < alt_ranges[1])
                add_vals = df_timestep[gnss_const][signal_family][metric][::plot_inv][alt_mask]
                plot_vals[metric] = np.concatenate([plot_vals.get(metric, np.array([])), add_vals])

        # Plotting
        for k, metric in enumerate(metrics_plot):
            # In axis 0, plot histogram
            bins_metric = np.linspace(range_metrics[k][0], range_metrics[k][1], 100)
            axes[k].hist(
                plot_vals[metric] * 1000,  # convert to mm
                bins=bins_metric,
                color=colors[i],
                alpha=0.5,
                density=True,
                label=f"{alt_ranges[0]}-{alt_ranges[1]} km",
            )
            axes[k].set_title(f"Histogram of {ylabels[k]}", fontsize=fs_title, fontweight="bold")
            axes[k].set_xlabel(ylabels[k], fontsize=fs)
            axes[k].set_ylabel("Density", fontsize=fs)
            axes[k].grid(True)
            axes[k].legend(loc="upper right", fontsize=fs_legend)
            axes[k].set_xlim(range_metrics[k])
            axes[k].tick_params(axis="both", labelsize=fs_ticks)

            is_not_outliers = np.abs(plot_vals[metric] * 1000) < 10  # 10 mm threshold
            plot_vals_without_outliers = plot_vals[metric][is_not_outliers]
            print(
                f"[TDCP {ylabels[k]}] Altitude Range: {alt_ranges[0]}-{alt_ranges[1]} km  mean: ",
                np.mean(plot_vals_without_outliers) * 1000,
                " mm,  std: ",
                np.std(plot_vals_without_outliers) * 1000,
                " mm",
            )

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300)
    plt.show()


def plot_iono_free_ure_errors(df_timestep, gnss_consts, plot_inv=10, filename=None):
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    fs = 16
    fs_legend = 14
    fs_ticks = 14
    fs_title = 16
    colors = ["blue", "orange", "green", "purple"]

    # Histogram of ionofree ure errors
    for i, gnss_const in enumerate(gnss_consts):
        ure_errors = df_timestep[gnss_const]["L1L5"]["ure_error"][::plot_inv]

        bins_ure = np.linspace(-3, 3, 100)  # remove outliers
        ax.hist(
            ure_errors, bins=bins_ure, color=colors[i], alpha=0.7, label=gnss_const, density=False
        )
        ax.set_title(
            "Histogram of L1L5 Ionofree Ephemeris Errors", fontsize=fs_title, fontweight="bold"
        )
        ax.set_xlabel("Ionofree Ephemeris Error (m)", fontsize=fs)
        ax.set_ylabel("Density", fontsize=fs)
        ax.grid(True)
        ax.legend(fontsize=fs_legend)

        is_not_outliers = np.abs(ure_errors) < 5  # 1 m threshold
        ure_errors_without_outliers = ure_errors[is_not_outliers]
        print(
            f"[L1L5 Ionofree URE] GNSS: {gnss_const}  mean: ",
            np.mean(ure_errors_without_outliers),
            " m,  std: ",
            np.std(ure_errors_without_outliers),
            " m",
        )

    if filename is not None:
        plt.savefig(filename, dpi=300)
    plt.show()


def plot_iono_delays_altitude(df_timestep, plot_inv=10, filename=None):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    fs = 16
    fs_legend = 14
    fs_ticks = 14
    fs_title = 16
    colors_gnss = ["blue", "orange", "green"]
    colors_freq = ["red", "purple"]
    signal_labels = ["L1/E1", "L5/E5a"]
    signal_families = [1, 5]

    for j, signal_idx in enumerate(signal_families):
        plot_data = {}
        for i, gnss_const in enumerate(df_timestep.keys()):
            # Axes 1: Min Alt vs Iono Delay Scatter Plot
            delays_iono = df_timestep[gnss_const][signal_idx]["pr_iono_error"][::plot_inv]
            min_alts = df_timestep[gnss_const][signal_idx]["min_alt_pr"][::plot_inv]  # in km
            tdcp_iono = df_timestep[gnss_const][signal_idx]["tdcp_iono_error"][::plot_inv]
            tdcp_min_alts = df_timestep[gnss_const][signal_idx]["min_alt_tdcp"][::plot_inv]  # in km

            plot_data["delays_iono"] = np.concatenate(
                [plot_data.get("delays_iono", np.array([])), delays_iono]
            )
            plot_data["min_alts"] = np.concatenate(
                [plot_data.get("min_alts", np.array([])), min_alts]
            )
            plot_data["tdcp_iono"] = np.concatenate(
                [plot_data.get("tdcp_iono", np.array([])), tdcp_iono]
            )
            plot_data["tdcp_min_alts"] = np.concatenate(
                [plot_data.get("tdcp_min_alts", np.array([])), tdcp_min_alts]
            )

        # Axes1: Min alt vs Delay scatter plot
        axes[0].scatter(
            plot_data["min_alts"],
            plot_data["delays_iono"],  # convert to m
            color=colors_freq[j],
            alpha=0.3,
            label=signal_labels[j],
        )
        axes[0].set_title(f"Pseudorange Ionospheric Delays", fontsize=fs_title, fontweight="bold")
        axes[0].set_xlabel("Tangential Altitude (km)", fontsize=fs)
        axes[0].set_ylabel("Ionospheric Delay (m)", fontsize=fs)
        axes[0].grid(True)
        axes[0].set_xlim([0, 20000])
        axes[0].set_yscale("log")
        axes[0].legend(fontsize=fs_legend)
        axes[0].tick_params(axis="both", labelsize=fs_ticks)

        # Axes 2: Min Alt vs TDCP Iono Delay Scatter Plot
        axes[1].scatter(
            plot_data["tdcp_min_alts"],
            plot_data["tdcp_iono"],  # convert to mm
            color=colors_freq[j],
            alpha=0.3,
            label=signal_labels[j],
        )
        axes[1].set_title(f"TDCP Ionospheric Delays", fontsize=fs_title, fontweight="bold")
        axes[1].set_xlabel("Tangential Altitude (km)", fontsize=fs)
        axes[1].set_ylabel("TDCP Ionospheric Delay (m)", fontsize=fs)
        axes[1].grid(True)
        axes[1].set_xlim([0, 20000])
        axes[1].set_yscale("log")
        axes[1].legend(fontsize=fs_legend)
        axes[1].tick_params(axis="both", labelsize=fs_ticks)

    for i, gnss_const in enumerate(df_timestep.keys()):
        # Axes 3: Min Alt vs Iono-free Delay Scatter Plot
        delays_iono_ionofree = df_timestep[gnss_const]["L1L5"]["delays_iono"]
        min_alts_ionofree = df_timestep[gnss_const]["L1L5"]["min_alt"]  # in km
        axes[2].scatter(
            min_alts_ionofree,
            delays_iono_ionofree,  # convert to m
            color=colors_gnss[i],
            alpha=0.5,
            label=gnss_const,
        )
        axes[2].set_title(f" L1L5 Ionofree", fontsize=fs_title, fontweight="bold")
        axes[2].set_xlabel("Tangential Altitude (km)", fontsize=fs)
        axes[2].set_ylabel("Ionospheric Delay (m)", fontsize=fs)
        axes[2].grid(True)
        axes[2].set_xlim([0, 2000])
        axes[2].set_yscale("log")
        axes[2].set_ylim([1e-7, 1e2])
        axes[2].set_yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2])
        axes[2].legend(fontsize=fs_legend)
        axes[2].tick_params(axis="both", labelsize=fs_ticks)

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300)
    plt.show()


def plot_noise_altitude(timestep_df, plot_inv=10, filename=None):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    fs = 16
    fs_legend = 14
    fs_title = 16
    fs_ticks = 14
    colors_gnss = ["blue", "orange", "green"]
    colors_freq = ["purple", "red"]
    signal_labels = ["L5/E5a", "L1/E1"]
    signal_families = [5, 1]

    for j, signal_idx in enumerate(signal_families):
        plot_data = {}
        for i, gnss_const in enumerate(timestep_df.keys()):
            # Axes 1: Min Alt vs PR Noise Scatter Plot
            pr_noise = timestep_df[gnss_const][signal_idx]["pr_noise"][::plot_inv]
            min_alts_pr = timestep_df[gnss_const][signal_idx]["min_alt_pr"][::plot_inv]  # in km
            tdcp_noise = timestep_df[gnss_const][signal_idx]["tdcp_noise"][::plot_inv]
            min_alts_tdcp = timestep_df[gnss_const][signal_idx]["min_alt_tdcp"][::plot_inv]  # in km

            plot_data["pr_noise"] = np.concatenate(
                [plot_data.get("pr_noise", np.array([])), pr_noise]
            )
            plot_data["min_alts_pr"] = np.concatenate(
                [plot_data.get("min_alts_pr", np.array([])), min_alts_pr]
            )
            plot_data["tdcp_noise"] = np.concatenate(
                [plot_data.get("tdcp_noise", np.array([])), tdcp_noise]
            )
            plot_data["min_alts_tdcp"] = np.concatenate(
                [plot_data.get("min_alts_tdcp", np.array([])), min_alts_tdcp]
            )

        # Axes1: Min alt vs PR Noise scatter plot
        axes[0].scatter(
            plot_data["min_alts_pr"],
            plot_data["pr_noise"],
            color=colors_freq[j],
            alpha=0.1,
            label=signal_labels[j],
        )
        axes[0].set_title(f"Pseudorange Receiver Noise", fontsize=fs_title, fontweight="bold")
        axes[0].set_xlabel("Tangential Altitude (km)", fontsize=fs)
        axes[0].set_ylabel("Pseudorange Noise (mm)", fontsize=fs)
        axes[0].grid(True)
        axes[0].set_xlim([0, 20000])
        axes[0].legend(fontsize=fs_legend)
        axes[0].tick_params(axis="both", labelsize=fs_ticks)

        # Axes 2: Min Alt vs TDCP Noise Scatter Plot
        axes[1].scatter(
            plot_data["min_alts_tdcp"],
            plot_data["tdcp_noise"] * 1000,  # convert to mm
            color=colors_freq[j],
            alpha=0.3,
            label=signal_labels[j],
        )
        axes[1].set_title(f"TDCP Receiver Noise", fontsize=fs_title, fontweight="bold")
        axes[1].set_xlabel("Tangential Altitude (km)", fontsize=fs)
        axes[1].set_ylabel("TDCP Noise (mm)", fontsize=fs)
        axes[1].grid(True)
        axes[1].set_xlim([0, 20000])
        axes[1].legend(fontsize=fs_legend)
        axes[1].tick_params(axis="both", labelsize=fs_ticks)

    # Iono-free noise plot
    for i, gnss_const in enumerate(timestep_df.keys()):
        # Axes 3: Min Alt vs Iono-free PR Noise Scatter Plot
        pr_noise_ionofree = timestep_df[gnss_const]["L1L5"]["pr_noise"][::plot_inv]
        min_alts_ionofree = timestep_df[gnss_const]["L1L5"]["min_alt"][::plot_inv]  #
        axes[2].scatter(
            min_alts_ionofree, pr_noise_ionofree, color=colors_gnss[i], alpha=0.5, label=gnss_const
        )
        axes[2].set_title(f"L1L5 Ionofree Receiver Noise", fontsize=fs_title, fontweight="bold")
        axes[2].set_xlabel("Tangential Altitude (km)", fontsize=fs)
        axes[2].set_ylabel("Iono-free Noise (m)", fontsize=fs)
        axes[2].grid(True)
        axes[2].set_xlim([0, 20000])
        axes[2].legend(fontsize=fs_legend)
        axes[2].tick_params(axis="both", labelsize=fs_ticks)

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300)
    plt.show()
