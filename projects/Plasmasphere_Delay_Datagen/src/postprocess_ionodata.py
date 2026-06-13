import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import pylupnt as pnt
from tqdm import tqdm
from src.clock_models import ClockNoise
import time
import pandas as pd
from multiprocessing import get_context, cpu_count
from scipy.interpolate import CubicSpline


def to_floatarray(str_array):
    # first split with spaces
    str_split = str_array.strip().split()
    if str_split[0] == "[":
        arr1 = float(str_split[1])
        arr2 = float(str_split[2])
        # if split[3] contains closing parenthesis, remove it
        if str_split[3].endswith("]"):
            arr3 = float(str_split[3][:-1])
        else:
            arr3 = float(str_split[3])
    else:
        arr1 = float(str_split[0][1:])  # remove the opening parenthesis
        arr2 = float(str_split[1])
        # if split[2] contains closing parenthesis, remove it
        if str_split[2].endswith("]"):
            arr3 = float(str_split[2][:-1])
        else:
            arr3 = float(str_split[2])

    float_list = [arr1, arr2, arr3]

    return np.array(float_list)


def plot_labels_data(df, gnss_consts_plot, inv=1, filename=None):
    min_alts = df["min_alt"].values[::inv]  # convert to km
    sigma_ranges = df["sigma_range"].values[::inv] * 1000  # to meters
    cn0s = df["cn0"].values[::inv]
    gnss_consts = df["gnss_const"].values[::inv]
    pos_tx = df["pos_tx"].values[::inv]
    vel_tx = df["vel_tx"].values[::inv]
    pos_rx = df["pos_rx"].values[::inv]
    clock_tx = df["clockbias_tx"].values[::inv]
    pos_tx_ephem = df["pos_tx_ephem"].values[::inv]
    vel_tx_ephem = df["vel_tx_ephem"].values[::inv]
    clock_tx_ephem = df["clockbias_ephem"].values[::inv]

    pos_tx_diff = (
        np.array(
            [
                np.linalg.norm(pos_tx[i] - pos_rx[i]) - np.linalg.norm(pos_tx_ephem[i] - pos_rx[i])
                for i in range(len(pos_tx))
            ]
        )
        * 1000
    )
    vel_tx_diff = (
        np.array([np.linalg.norm(vel_tx[i] - vel_tx_ephem[i]) for i in range(len(vel_tx))]) * 1e6
    )

    C_ms = 299792458.0  # speed of light in m/s
    clock_tx_diff = (clock_tx - clock_tx_ephem) * C_ms
    ure_tx_diff = pos_tx_diff + clock_tx_diff
    print(clock_tx_diff)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    colors = ["blue", "orange", "purple", "cyan", "magenta", "yellow"]
    fs = 14

    # set tick param size
    for ax in axes:
        ax.tick_params(axis="both", labelsize=12)  # 'x', 'y', or 'both'  # font size

    for i, gnss_const in enumerate(gnss_consts_plot):
        mask = gnss_consts == gnss_const
        axes[0].hist(min_alts[mask], bins=30, color=colors[i], alpha=0.7, label=gnss_const)
        axes[0].set_title("Histogram of Tangential Altitude (min_alt)")
        axes[0].set_xlabel("Altitude (km)", fontsize=fs)
        axes[0].set_ylabel("Frequency", fontsize=fs)
        axes[0].grid(True)
        axes[0].legend()
        axes[1].hist(sigma_ranges[mask], bins=30, color=colors[i], alpha=0.7, label=gnss_const)
        axes[1].set_title("Histogram of DLL error (sigma_range)")
        axes[1].set_xlabel("Sigma Range (m)", fontsize=fs)
        axes[1].set_ylabel("Frequency", fontsize=fs)
        axes[1].grid(True)
        axes[1].legend()
        axes[2].hist(cn0s[mask], bins=30, color=colors[i], alpha=0.7, label=gnss_const)
        axes[2].set_title("Histogram of cn0")
        axes[2].set_xlabel("cn0 (dB-Hz)", fontsize=fs)
        axes[2].set_ylabel("Frequency", fontsize=fs)
        axes[2].grid(True)
        axes[2].legend()

        bins_pos = np.linspace(-5, 5, 100)  # remove outliers
        axes[3].hist(pos_tx_diff[mask], bins=bins_pos, color=colors[i], alpha=0.7, label=gnss_const)
        axes[3].set_title("Histogram of Transmitter Position Difference")
        axes[3].set_xlabel("Position Difference (m)", fontsize=fs)
        axes[3].set_ylabel("Frequency", fontsize=fs)
        axes[3].grid(True)
        axes[3].legend()

        bins_clk = np.linspace(-5, 5, 100)  # remove outliers
        axes[4].hist(
            clock_tx_diff[mask], bins=bins_clk, color=colors[i], alpha=0.7, label=gnss_const
        )
        axes[4].set_title("Histogram of Transmitter Clock Bias Difference")
        axes[4].set_xlabel("Clock Bias Difference (m)", fontsize=fs)
        axes[4].set_ylabel("Frequency", fontsize=fs)
        axes[4].grid(True)
        axes[4].legend()

        bins_vel = np.linspace(-5, 5, 100)  # remove outliers
        axes[5].hist(ure_tx_diff[mask], bins=bins_vel, color=colors[i], alpha=0.7, label=gnss_const)
        axes[5].set_title("Histogram of Transmitter URE Difference")
        axes[5].set_xlabel("URE Difference (m)", fontsize=fs)
        axes[5].set_ylabel("Frequency", fontsize=fs)
        axes[5].grid(True)
        axes[5].legend()

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename)
    plt.show()


def plot_delay_data(df, gnss_consts_plot, freqs_plot=[1, 5]):
    total_delay_m = df["total_delay_m"].values
    tec_delay_m = df["tec_delay_m"].values
    min_alts = df["min_alt"].values  # in km
    gnss_consts = df["gnss_const"].values
    sigma_ranges = df["sigma_range"].values * 1000  # to meters
    tec = df["tecu"].values
    dist_bend_m = df["dist_bend_m"].values

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    colors = ["blue", "orange", "purple", "cyan", "magenta", "yellow"]

    axes = axes.flatten()

    i = 0

    fs = 14
    for ax in axes:
        ax.tick_params(axis="both", labelsize=12)  # 'x', 'y', or 'both'  # font size

    for ii, gnss_const in enumerate(gnss_consts_plot):
        for jj, freq in enumerate(freqs_plot):
            if gnss_const == "GPS" and freq == 1:
                freq_label = "L1"
                freq_label_plot = "L1"
            elif gnss_const == "GPS" and freq == 5:
                freq_label = "L5"
                freq_label_plot = "L5"
            elif gnss_const == "GALILEO" and freq == 1:
                freq_label = "E1"
                freq_label_plot = "E1"
            elif gnss_const == "GALILEO" and freq == 5:
                freq_label = "E5a"
                freq_label_plot = "E5a"
            elif gnss_const == "QZSS" and freq == 1:
                freq_label = "L1"
                freq_label_plot = "L1"
            elif gnss_const == "QZSS" and freq == 5:
                freq_label = "L5"
                freq_label_plot = "L5"

            mask = (gnss_consts == gnss_const) & (df["signal"] == freq_label)
            label_prefix = f"{gnss_const} {freq_label_plot}"
            axes[0].hist(tec_delay_m[mask], bins=30, color=colors[i], alpha=0.4, label=label_prefix)
            axes[0].set_title("Histogram of TEC Delay")
            axes[0].set_xlabel("TEC Delay (m)", fontsize=fs)
            axes[0].set_ylabel("Frequency", fontsize=fs)
            axes[0].grid(True)
            axes[0].set_yscale("log")
            axes[0].legend()

            # minalt vs total delay
            axes[1].scatter(
                min_alts[mask], total_delay_m[mask], color=colors[i], alpha=0.3, label=label_prefix
            )
            axes[1].set_title("Tagential Altitude vs Total Delay")
            axes[1].set_xlabel("Tagential Altitude (km)", fontsize=fs)
            axes[1].set_ylabel("Total Delay (m)", fontsize=fs)
            axes[1].grid(True)
            axes[1].set_yscale("log")
            axes[1].set_xlim(0, 20000)
            axes[1].legend()

            # range sigma vs total delay
            axes[2].scatter(
                sigma_ranges[mask],
                total_delay_m[mask],
                color=colors[i],
                alpha=0.3,
                label=label_prefix,
            )
            axes[2].set_title("range sigma vs Total Delay")
            axes[2].set_xlabel("range sigma (m)", fontsize=fs)
            axes[2].set_ylabel("Total Delay (m)", fontsize=fs)
            axes[2].grid(True)
            axes[2].set_yscale("log")
            axes[2].legend()

            # minalt vs range sigma
            axes[3].scatter(
                min_alts[mask], sigma_ranges[mask], color=colors[i], alpha=0.3, label=label_prefix
            )
            axes[3].set_title("Tagential Altitude vs range sigma")
            axes[3].set_xlabel("Tagential Altitude (km)", fontsize=fs)
            axes[3].set_ylabel("range sigma (m)", fontsize=fs)
            axes[3].grid(True)
            axes[3].legend()
            axes[3].set_xlim(0, 20000)

            # minalt vs tec
            axes[4].scatter(
                min_alts[mask], tec_delay_m[mask], color=colors[i], alpha=0.3, label=label_prefix
            )
            axes[4].set_title("Tagential Altitude vs tec delay")
            axes[4].set_xlabel("Tagential Altitude (km)", fontsize=fs)
            axes[4].set_ylabel("TEC Delay (m)", fontsize=fs)
            axes[4].grid(True)
            axes[4].set_yscale("log")
            axes[4].legend()
            axes[4].set_xlim(0, 20000)

            # minalt vs dist_bend_m
            axes[5].scatter(
                min_alts[mask], dist_bend_m[mask], color=colors[i], alpha=0.3, label=label_prefix
            )
            axes[5].set_title("Tagential Altitude vs Distance Bend")
            axes[5].set_xlabel("Tagential Altitude (km)", fontsize=fs)
            axes[5].set_ylabel("Distance Bend (m)", fontsize=fs)
            axes[5].grid(True)
            axes[5].set_yscale("log")
            axes[5].legend()
            axes[5].set_xlim(0, 2000)

            i = i + 1

    plt.tight_layout()
    plt.show()


def load_raytrace_csv_data(
    sim_params, datadir, raytracedir, savedir, overwrite=False, load_correction=True
):
    ymdh = sim_params["epoch_ymdh"]
    lcrns_idx = sim_params["lcrns_idx"]
    signal_family = sim_params["signal_family"]
    rz12 = sim_params["rz12"]
    gnss_const = sim_params["gnss_const"]
    kp = sim_params.get("kp", None)
    n_orbit = sim_params.get("n_orbit", 3)
    dt = sim_params.get("dt", 10)
    dtrt = sim_params.get("dt_raytrace", 60)

    orbit_dt_str = "norbit_{}_dt_{}s_dtrt_{}s".format(int(n_orbit), int(dt), int(dtrt))

    epoch_dict = {
        "year": ymdh[0],
        "month": ymdh[1],
        "day": ymdh[2],
        "hour": ymdh[3],
        "minute": 0,
        "second": 0,
    }
    epoch_str = "{year}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}_{second:02d}".format(
        **epoch_dict
    )

    outfilename = "ionodata_raytrace_{}_sat_{}_{}_signal_{}_rz12_{:.1f}_kp_{:.1f}.pkl".format(
        epoch_str, lcrns_idx, gnss_const, signal_family, rz12, kp
    )
    savedir = os.path.join(savedir, orbit_dt_str)
    updated_pickle_filename = os.path.join(savedir, outfilename)

    if not os.path.exists(savedir):
        os.makedirs(savedir)

    if os.path.exists(updated_pickle_filename) and not overwrite:
        print(f"Updated CSV file {outfilename} already exists in {savedir}. Loading it.")
        df_existing = pd.read_pickle(updated_pickle_filename)
        return df_existing

    # load the pickle file with full information -------------------------------------------------
    pkl_labels_filename = os.path.join(
        datadir,
        orbit_dt_str,
        "ionodata_full_{}_sat_{}_{}_signal_{}.pkl".format(
            epoch_str, lcrns_idx, gnss_const, signal_family
        ),
    )

    if not os.path.exists(pkl_labels_filename):
        raise FileNotFoundError(f"Pickle file {pkl_labels_filename} does not exist.")

    # Load the initial ionodata CSV
    print(f"Loading pickle file {pkl_labels_filename}...")
    df = pd.read_pickle(pkl_labels_filename)
    print(f"Pickle file {pkl_labels_filename} loaded successfully.")
    print("Number of rows in the dataframe:", len(df))

    # load the shortened raytrace pickle file -------------------------------------------------
    pkl_labels_filename_rt = os.path.join(
        datadir,
        orbit_dt_str,
        "ionodata_short_{}_sat_{}_{}_signal_{}.pkl".format(
            epoch_str, lcrns_idx, gnss_const, signal_family
        ),
    )
    if not os.path.exists(pkl_labels_filename_rt):
        raise FileNotFoundError(
            f"Short raytrace pickle file {pkl_labels_filename_rt} does not exist."
        )
    print(f"Loading short raytrace pickle file {pkl_labels_filename_rt}...")
    df_rt = pd.read_pickle(pkl_labels_filename_rt)
    print(f"Short raytrace pickle file {pkl_labels_filename_rt} loaded successfully.")
    print("Number of rows in the shortened raytrace dataframe:", len(df_rt))
    row_full = df_rt["row_full"].values

    # overwrite the delay datas (always load signa 1 case) -------------------------------------------------
    orbit_dtrt_str = "norbit_{}_dt_{}s_dtrt_{}s".format(int(n_orbit), int(dt), int(dtrt))
    ray_trace_csv_dir = os.path.join(
        raytracedir,
        orbit_dtrt_str,
        "raytrace_results_sat_{}_{}_signal_{}_rz12_{:.0f}_kp_{:.0f}".format(
            lcrns_idx, gnss_const, signal_family, rz12, kp
        ),
    )

    if not os.path.exists(ray_trace_csv_dir):
        raise FileNotFoundError(f"Ray trace CSV directory {ray_trace_csv_dir} does not exist.")
    else:
        # Load and overwrite delay data from ray trace CSV files
        # list all csv files in the directory
        ray_trace_files = [f for f in os.listdir(ray_trace_csv_dir) if f.endswith(".csv")]
        print("Number of ray trace CSV files found:", len(ray_trace_files))
        raytrace_cols = [
            "total_delay_m",
            "tecu",
            "tec_delay_m",
            "second_delay_m",
            "third_delay_m",
            "dist_bend_m",
            "tec_delay_bend_m",
            "max_sep_line_m",
            "final_pos_err_m",
        ]
        for file in tqdm(ray_trace_files, desc="Loading ray trace CSV files"):
            # extract row number from filename (raytrace_result_{row}.csv)
            row_rt = int(file.split("_")[-1].split(".")[0])

            if row_rt >= len(row_full):
                print(
                    f"Warning: row raytrace {row_rt} exceeds row_full length {len(row_full)}. Skipping."
                )
                continue
            row_num = row_full[row_rt]
            if row_num >= len(df):
                print(f"Warning: row_num {row_num} exceeds dataframe length {len(df)}. Skipping.")
                continue
            df_ray = pd.read_csv(os.path.join(ray_trace_csv_dir, file))
            # replace columns
            # total_delay_m,tecu,tec_delay_m,second_delay_m,third_delay_m,dist_bend_m,tec_delay_bend_m,max_sep_line_m,final_pos_err_m
            data = df_ray[raytrace_cols].iloc[0]
            # replace the corresponding row in df
            for col in raytrace_cols:
                df.at[row_num, col] = data[col]

        # Show the rows that are missing
        # if len(row_nums) < len(df):
        #     missing_rows = set(range(len(df))) - set(row_nums)
        #     print(f"Warning: Missing ray trace data for rows: {sorted(missing_rows)}")
        #     print

    # Load and overwrite delay data with corrections -------------------------------------------------
    correct_csv_dir = os.path.join(
        raytracedir,
        orbit_dtrt_str,
        "raytrace_results_correct_sat_{}_{}_signal_{}_rz12_{:.0f}_kp_{:.0f}".format(
            lcrns_idx, gnss_const, signal_family, rz12, kp
        ),
    )

    if not os.path.exists(correct_csv_dir):
        print(f"Corrected ray trace CSV directory {correct_csv_dir} does not exist.")
    elif not load_correction:
        print(
            f"Skipping loading corrected ray trace CSV files from {correct_csv_dir} as load_correction is False."
        )
    else:
        print(f"Loading corrected ray trace CSV files from {correct_csv_dir}.")
        correct_files = [f for f in os.listdir(correct_csv_dir) if f.endswith(".csv")]
        for file in tqdm(correct_files, desc="Loading corrected ray trace CSV files"):
            # extract row number from filename (raytrace_result_{row}.csv)
            row_rt = int(file.split("_")[-1].split(".")[0])

            if row_rt >= len(row_full):
                print(
                    f"Warning: row raytrace {row_rt} exceeds row_full length {len(row_full)}. Skipping."
                )
                continue
            row_num = row_full[row_rt]
            if row_num >= len(df):
                print(f"Warning: row_num {row_num} exceeds dataframe length {len(df)}. Skipping.")
                continue
            df_correct = pd.read_csv(os.path.join(correct_csv_dir, file))
            # replace columns
            data = df_correct[raytrace_cols].iloc[0]
            # replace the corresponding row in df
            for col in raytrace_cols:
                df.at[row_num, col] = data[col]

    # Interpolate missing raytrace data -------------------------------------------------
    print("Interpolating missing raytrace data...")
    df = interpolate_missing_raytrace_data(df, dt, dtrt, raytrace_cols)

    # Save the updated dataframe to a new CSV file ------------------------------------------
    # to pickle
    print(f"Saving updated dataframe to pickle file {updated_pickle_filename}...")
    df.to_pickle(updated_pickle_filename)
    print(f"Updated pickle file saved to {updated_pickle_filename}.")

    print("Conversion to CSV started...")
    updated_csv_filename = updated_pickle_filename.replace(".pkl", ".csv")
    df.to_csv(updated_csv_filename, index=False)
    print(f"Updated CSV file saved to {updated_csv_filename}.")

    return df


def interpolate_missing_raytrace_data(df, dt, dtrt, raytrace_cols):
    """
    Interpolate missing raytrace data in the dataframe.

    :param df: DataFrame containing raytrace data
    :param dt: Time difference used for interpolation
    :param raytrace_cols: List of columns to interpolate
    """

    df_new = df.copy()
    missing_rows = df[df["total_delay_m"] == 0].index.tolist()
    prns = df["prn"].values
    unique_prns = np.unique(prns)
    print("Unique PRNs:", unique_prns)

    rt_start_idx = 1  # raytrace start time index
    rt_inv = int(dtrt / dt)

    def prev_rt_tidx(tidx):
        return rt_inv * np.floor((tidx - rt_start_idx) / rt_inv) + rt_start_idx

    def next_rt_tidx(tidx):
        return rt_inv * np.ceil((tidx - rt_start_idx) / rt_inv) + rt_start_idx

    if len(missing_rows) > 0:
        for prn in tqdm(unique_prns, desc="Interpolating PRNs"):
            mask_prn = prns == prn
            df_prn = df[mask_prn]
            is_raytrace = df_prn["total_delay_m"] != 0
            df_prn_missing = df_prn[~is_raytrace]
            df_prn_not_missing = df_prn[is_raytrace]

            # extract tspan and get points where prn is also tracked at before and after
            tidx_prn_missing = df_prn_missing["tidx"].values
            tidx_prn_not_missing = df_prn_not_missing["tidx"].values

            # first, find missing rows that have valid data before and after
            if len(df_prn_not_missing) < 2:
                continue

            # For each row with missing raytrace data, find the previous and next raytrace indices
            tidx_prn_missing_next_rtidx = next_rt_tidx(tidx_prn_missing)
            tidx_prn_prev_rtidx = prev_rt_tidx(tidx_prn_missing)

            # Check which missing rows have previous and next raytrace indices available
            prn_prev_rt_in = np.isin(tidx_prn_prev_rtidx, tidx_prn_not_missing)
            prn_next_rt_in = np.isin(tidx_prn_missing_next_rtidx, tidx_prn_not_missing)

            # Case 1: middle values (both before and after exist)
            df_prn_middle = df_prn_missing[prn_prev_rt_in & prn_next_rt_in]

            # Case 2: starting values (only after exist)
            df_prn_start = df_prn_missing[~prn_prev_rt_in & prn_next_rt_in]

            # Case 3: ending values (only before exist)
            df_prn_end = df_prn_missing[prn_prev_rt_in & ~prn_next_rt_in]

            # Case 4: isolated values (neither before nor after exist)
            df_prn_isolated = df_prn_missing[~prn_prev_rt_in & ~prn_next_rt_in]

            print(
                "PRN:",
                prn,
                " | Missing total: ",
                len(df_prn_missing),
                "  Middle: ",
                len(df_prn_middle),
                "  Starting: ",
                len(df_prn_start),
                "  Ending: ",
                len(df_prn_end),
                "  Isolated:",
                len(df_prn_isolated),
            )

            if len(df_prn_middle) == 0:
                continue

            # x: tspan where this PRN has valid raytrace
            x = pd.to_numeric(df_prn_not_missing["tspan"], errors="coerce").to_numpy(dtype=float)
            # x_new: tspan of missing rows
            x_new = pd.to_numeric(df_prn_middle["tspan"], errors="coerce").to_numpy(dtype=float)

            for col in raytrace_cols:
                # y: corresponding column values
                y = pd.to_numeric(df_prn_not_missing[col], errors="coerce").to_numpy(dtype=float)

                # optional: drop NaNs from x/y before interp
                valid = ~np.isnan(x) & ~np.isnan(y)
                if valid.sum() < 2:  # not enough points to interpolate
                    continue

                # use_cubic = False  # To use cubic, we need to run for each seperate blocks
                # if use_cubic and valid.sum() >= 4:
                #     # use cubic interpolation
                #     cs = CubicSpline(x[valid], y[valid])
                #     df_new.loc[df_prn_middle.index, col] = cs(x_new)
                # else:
                #     # Case 1: middle values -> do linear interpolation
                #     df_new.loc[df_prn_middle.index, col] = np.interp(x_new, x[valid], y[valid])
                df_new.loc[df_prn_middle.index, col] = np.interp(x_new, x[valid], y[valid])

                # Case 2: starting values -> extrapolate
                if len(df_prn_start) > 0:
                    for idx in df_prn_start.index:
                        tidx_next = next_rt_tidx(df_prn_start.at[idx, "tidx"])
                        next_value = df_prn_not_missing.loc[
                            df_prn_not_missing["tidx"] == tidx_next, col
                        ].values

                        tidx_next2 = tidx_next + rt_inv
                        next_value2 = df_prn_not_missing.loc[
                            df_prn_not_missing["tidx"] == tidx_next2, col
                        ].values

                        if len(next_value) > 0 and len(next_value2) > 0:
                            tspan_val = df_prn_start.at[idx, "tspan"]
                            tspan_next = df_prn_not_missing.loc[
                                df_prn_not_missing["tidx"] == tidx_next, "tspan"
                            ].values
                            tspan_next2 = df_prn_not_missing.loc[
                                df_prn_not_missing["tidx"] == tidx_next2, "tspan"
                            ].values

                            val = np.interp(
                                tspan_val,
                                [tspan_next[0], tspan_next2[0]],
                                [next_value[0], next_value2[0]],
                            )
                            df_new.at[idx, col] = val
                        elif len(next_value) > 0:
                            df_new.at[idx, col] = next_value[0]

                # Case 3: ending values -> use previous value
                if len(df_prn_end) > 0:
                    for idx in df_prn_end.index:
                        tidx_prev = prev_rt_tidx(df_prn_end.at[idx, "tidx"])
                        prev_value = df_prn_not_missing.loc[
                            df_prn_not_missing["tidx"] == tidx_prev, col
                        ].values

                        tidx_prev2 = tidx_prev - rt_inv
                        prev_value2 = df_prn_not_missing.loc[
                            df_prn_not_missing["tidx"] == tidx_prev2, col
                        ].values
                        tspan_prev = df_prn_not_missing.loc[
                            df_prn_not_missing["tidx"] == tidx_prev, "tspan"
                        ].values
                        tspan_prev2 = df_prn_not_missing.loc[
                            df_prn_not_missing["tidx"] == tidx_prev2, "tspan"
                        ].values

                        if len(prev_value) > 0 and len(prev_value2) > 0:
                            val = np.interp(
                                df_prn_end.at[idx, "tspan"],
                                [tspan_prev[0], tspan_prev2[0]],
                                [prev_value[0], prev_value2[0]],
                            )
                            df_new.at[idx, col] = val
                        elif len(prev_value) > 0:
                            df_new.at[idx, col] = prev_value[0]

    return df_new


def load_L1_L5_labels(sim_params, datadir_labels, datadir_raytrace, overwrite=False):
    ymdh = sim_params["epoch_ymdh"]
    lcrns_idx = sim_params["lcrns_idx"]
    rz12 = sim_params["rz12"]
    gnss_const = sim_params["gnss_const"]
    kp = sim_params["kp"]

    epoch_dict = {
        "year": ymdh[0],
        "month": ymdh[1],
        "day": ymdh[2],
        "hour": ymdh[3],
        "minute": 0,
        "second": 0,
    }
    epoch_str = "{year}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}_{second:02d}".format(
        **epoch_dict
    )

    # load labels Pickle for L1/E1 signal ---------------------------------------------------------------
    pickle_labels_filename_L1 = os.path.join(
        datadir_labels,
        "ionodata_short_{}_sat_{}_{}_signal_1.pkl".format(epoch_str, lcrns_idx, gnss_const),
    )

    if not os.path.exists(pickle_labels_filename_L1):
        raise FileNotFoundError(f"Pickle file {pickle_labels_filename_L1} does not exist.")

    # Load the initial ionodata Pickle
    df_L1 = pd.read_pickle(pickle_labels_filename_L1)
    print(f"Pickle file {pickle_labels_filename_L1} loaded successfully.")
    print("Number of rows in the L1 dataframe:", len(df_L1))

    # Load the labels CSV for L5/E5a signal ---------------------------------------------------------------
    pickle_labels_filename_L5 = os.path.join(
        datadir_labels,
        "ionodata_short_{}_sat_{}_{}_signal_5.pkl".format(epoch_str, lcrns_idx, gnss_const),
    )
    if not os.path.exists(pickle_labels_filename_L5):
        raise FileNotFoundError(f"Pickle file {pickle_labels_filename_L5} does not exist.")

    df_L5 = pd.read_pickle(pickle_labels_filename_L5)
    print(f"Pickle file {pickle_labels_filename_L5} loaded successfully.")
    print("Number of rows in the L5 dataframe:", len(df_L5))

    # Load the delay data CSV for L1 signal ---------------------------------------------------------------
    ray_trace_csv_dir_L1 = os.path.join(
        datadir_raytrace,
        "raytrace_results_sat_{}_{}_signal_1_rz12_{:.0f}_kp_{:.0f}".format(
            lcrns_idx, gnss_const, rz12, kp
        ),
    )

    if not os.path.exists(ray_trace_csv_dir_L1):
        raise FileNotFoundError(f"Ray trace CSV directory {ray_trace_csv_dir_L1} does not exist.")
    else:
        # Load and overwrite delay data from ray trace CSV files
        # list all csv files in the directory
        ray_trace_files = [f for f in os.listdir(ray_trace_csv_dir_L1) if f.endswith(".csv")]
        print("Number of ray trace CSV files found:", len(ray_trace_files))
        raytrace_cols = [
            "total_delay_m",
            "tecu",
            "tec_delay_m",
            "second_delay_m",
            "third_delay_m",
            "dist_bend_m",
            "tec_delay_bend_m",
            "max_sep_line_m",
            "final_pos_err_m",
        ]
        for file in tqdm(ray_trace_files, desc="Loading ray trace CSV files"):
            # extract row number from filename (raytrace_result_{row}.csv)
            row_num = int(file.split("_")[-1].split(".")[0])
            df_ray = pd.read_csv(os.path.join(ray_trace_csv_dir_L1, file))
            # replace columns
            # total_delay_m,tecu,tec_delay_m,second_delay_m,third_delay_m,dist_bend_m,tec_delay_bend_m,max_sep_line_m,final_pos_err_m
            data = df_ray[raytrace_cols].iloc[0]
            # replace the corresponding row in df
            for col in raytrace_cols:
                df_L1.at[row_num, col] = data[col]

    return df_L1, df_L5


def convert_L1_to_L5(sim_params, df_L1, df_L5, datadir_raytrace, overwrite=False):
    """Convert L1/E1 ephemeris data to L5/E5a ephemeris data using scaling factors.
    Args:
        sim_params (dict): Simulation parameters containing 'gnss_const', 'lcrns_idx', 'epoch_ymdh', and 'rz12'.
        datadir_raytrace (str): Directory path where the ray trace CSV files are stored.
        overwrite (bool): Whether to overwrite existing files. Default is False.
    Returns:
        pd.DataFrame: DataFrame containing the converted L5/E5a ephemeris data.
    """
    # Setup --------------------------------------------------------------------------------
    ymdh = sim_params["epoch_ymdh"]
    epoch_dict = {
        "year": ymdh[0],
        "month": ymdh[1],
        "day": ymdh[2],
        "hour": ymdh[3],
        "minute": 0,
        "second": 0,
    }
    epoch_str = "{year}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}_{second:02d}".format(
        **epoch_dict
    )

    # For each L5/E5a label, find the corresponding L1/E1 label and copy the ephemeris data ----------------
    match_num = 0
    L1_freq = 1575.42  # L1/E1 frequency in MHz
    L5_freq = 1176.45  # L5/E5a frequency in MHz
    L1_L5_fr2 = (L1_freq / L5_freq) ** 2
    L1_L5_fr3 = (L1_freq / L5_freq) ** 3
    L1_L5_fr4 = (L1_freq / L5_freq) ** 4

    datadir_save = os.path.join(
        datadir_raytrace,
        "raytrace_results_sat_{}_{}_signal_5_rz12_{:.0f}_kp_{:.0f}".format(
            sim_params["lcrns_idx"],
            sim_params["gnss_const"],
            sim_params["rz12"],
            sim_params["kp"],
        ),
    )

    os.makedirs(datadir_save, exist_ok=True)

    # Iterate over each row in df_L5 and find the corresponding L1/E1 data
    L1_tidx_row = 0

    for idx, row in tqdm(
        df_L5.iterrows(), total=len(df_L5), desc="Merging L5/E5a data with L1/E1 data"
    ):
        tidx = row["tidx"]

        # Update L1_tidx_row_start and L1_tidx_row_end to narrow down the search rang
        L1_rows_tidx = []
        for L1_tidx_row in range(L1_tidx_row, len(df_L1)):
            if df_L1.at[L1_tidx_row, "tidx"] == tidx:
                L1_rows_tidx.append(L1_tidx_row)
            if L1_tidx_row >= len(df_L1):
                break
            if df_L1.at[L1_tidx_row, "tidx"] > tidx:
                break

        # Find the matching row in df_L1 based on min_alt
        match_L1 = df_L1.loc[L1_rows_tidx]
        match_L1 = match_L1[(match_L1["min_alt"] == row["min_alt"])]

        if not match_L1.empty:
            # There should be only one matching row
            df_data = match_L1.iloc[0]
            # Copy ephemeris data from match_row to the current row in df
            df_L5_row = pd.DataFrame(
                {
                    "total_delay_m": 0.0,
                    "tecu": df_data["tecu"],
                    "tec_delay_m": df_data["tec_delay_m"] * L1_L5_fr2,
                    "second_delay_m": df_data["second_delay_m"] * L1_L5_fr3,
                    "third_delay_m": df_data["third_delay_m"] * L1_L5_fr4,
                    "dist_bend_m": 0.0,
                    "tec_delay_bend_m": 0,
                    "max_sep_line_m": 0,
                    "final_pos_err_m": 0,
                },
                index=[0],
            )
            df_L5_row["total_delay_m"] = (
                df_L5_row["tec_delay_m"] + df_L5_row["second_delay_m"] + df_L5_row["third_delay_m"]
            )
            match_num += 1

            # save the updated df_L5 to CSV
            df_L5_filename = os.path.join(datadir_save, "raytrace_result_{}.csv".format(idx))

            # save the updated df_L5 to CSV
            if not os.path.exists(df_L5_filename) or overwrite:
                # save only if file does not exist
                df_L5_row.to_csv(df_L5_filename, index=False)

    # print how many L5/E5a rows were matched with L1/E1 data
    print(f"Number of matched L5/E5a rows with L1/E1 data: {match_num}/{len(df_L5)}")


def generate_meas_timestep_df(sim_params, datapath, savepath, worker_id=0):
    """
    Generate Measurement Timestep DataFrame by aggregating raytrace data.
    """
    signals = sim_params["signal_families"]
    gnss_consts = sim_params["gnss_consts"]
    lcrns_idx = sim_params["lcrns_idx"]
    epoch_str = sim_params["epoch_str"]
    rz12 = sim_params["rz12"]
    kp = sim_params.get("kp", None)
    dt = sim_params.get("dt", 1.0)
    n_orbit = sim_params.get("n_orbit", 3)
    dt_rt = sim_params.get("dt_raytrace", 60)

    orbit_dt_str = "norbit_{}_dt_{}s_dtrt_{}s".format(int(n_orbit), int(dt), int(dt_rt))
    datapath = os.path.join(datapath, orbit_dt_str)
    savepath = os.path.join(savepath, orbit_dt_str)

    if not os.path.exists(savepath):
        os.makedirs(savepath)

    df_cols = [
        "tidx",
        "tspan",
        "t_tai",
        "acq_prns",
        "tracked_prns",
        "lost_prns",
        "num_tracked",
        "min_alts_km",
        "cn0",
        "sigma_range_m",
        "sigma_rangerate_ms",
        "sigma_carrier_m",
        "pos_tx_m",
        "vel_tx_ms",
        "clockbias_tx_s",
        "pos_rx_m",
        "vel_rx_ms",
        "pos_tx_ephem_m",
        "vel_tx_ephem_ms",
        "clockbias_ephem_s",
        "total_delay_m",
        "total_carrier_delay_m",
        "tec_delay_m",
        "second_delay_m",
        "third_delay_m",
        "dist_bend_m",
        "tec_delay_bend_m",
    ]

    df_outs = []

    # load all combinations of signals and gnss constellations
    for signal_family in signals:
        for gnss_const in gnss_consts:
            sim_params = {
                "gnss_const": gnss_const,
                "signal_family": signal_family,
                "lcrns_idx": lcrns_idx,
                "epoch_ymdh": [2025, 3, 1, 12],
                "rz12": rz12,
            }
            # Load pickle data
            filename = "ionodata_raytrace_{}_sat_{}_{}_signal_{}_rz12_{:.1f}_kp_{:.1f}.pkl".format(
                epoch_str, lcrns_idx, gnss_const, signal_family, rz12, kp
            )
            print("Loading:", filename)
            pickle_filename = os.path.join(datapath, filename)
            df = pd.read_pickle(pickle_filename)

            # replace to numpy array columns
            # array_cols = ["pos_tx", "vel_tx", "pos_rx", "vel_rx", "pos_tx_ephem", "vel_tx_ephem"]
            # for col in array_cols:
            #     df[col] = df[col].apply(to_floatarray)

            # construct final dataframe
            lent = int((df["tidx"].values.max() - df["tidx"].values.min())) + 1
            t_idxs_df = np.arange(df["tidx"].values.min(), df["tidx"].values.min() + lent)
            tspans_df = np.linspace(
                df["tspan"].values.min(), df["tspan"].values.min() + (lent - 1) * dt, lent
            )
            t_tais_df = np.linspace(
                df["t_tai"].values.min(), df["t_tai"].values.min() + (lent - 1) * dt, lent
            )

            tracked_prns = []
            df_out = pd.DataFrame(columns=df_cols, index=np.arange(lent))

            start_time = time.time()
            tidxs = df["tidx"].values
            start_idx = 0

            # load data for each time index --------------------------------------------------------
            print(f"[worker {worker_id}] Generating measurement timestep dataframe...")
            for ii, tidx in enumerate(t_idxs_df):

                if ii % 1000 == 0:
                    elapsed_time = time.time() - start_time
                    elapsed_hr = elapsed_time / 3600
                    elapsed_min = (elapsed_time % 3600) / 60
                    elapsed_sec = elapsed_time % 60
                    eta_sec = (lent - ii) * (elapsed_time / (ii + 1))
                    eta_hr = eta_sec / 3600
                    eta_min = (eta_sec % 3600) / 60
                    eta_sec_remain = eta_sec % 60
                    print(
                        f"[worker {worker_id}] Processing tidx {tidx} ({ii+1}/{lent}) - Elapsed: {int(elapsed_hr)}h {int(elapsed_min)}m {int(elapsed_sec)}s - ETA: {int(eta_hr)}h {int(eta_min)}m {int(eta_sec_remain)}s"
                    )

                # compute the row number
                # ii = tidx - df["tidx"].values[0] # zero-based index
                idxs = []
                while tidxs[start_idx] == tidx:
                    idxs.append(start_idx)
                    start_idx += 1
                    if start_idx >= len(tidxs):
                        break
                    if tidxs[start_idx] > tidx:
                        break

                df_tidx = df.iloc[idxs]
                df_out.at[ii, "tidx"] = tidx

                if len(df_tidx) == 0:
                    # print(f"Warning: No data for tidx {tidx}. Filling with empty values.")
                    df_out.at[ii, "tspan"] = tspans_df[ii]
                    df_out.at[ii, "t_tai"] = t_tais_df[ii]

                    for col in df_cols:
                        if (col != "tidx") and (col != "tspan") and (col != "t_tai"):
                            df_out.at[ii, col] = np.nan
                    continue

                tspan = df_tidx["tspan"].values[0]
                t_tai = df_tidx["t_tai"].values[0]
                df_out.at[ii, "tspan"] = tspan
                df_out.at[ii, "t_tai"] = t_tai
                # tracked PRNs
                curr_tracked_prns = df_tidx["prn"].tolist()
                df_out.at[ii, "tracked_prns"] = curr_tracked_prns
                # newly acquired PRNs
                acq_prns = [prn for prn in curr_tracked_prns if prn not in tracked_prns]
                df_out.at[ii, "acq_prns"] = acq_prns
                # lost PRNs
                lost_prns = [prn for prn in tracked_prns if prn not in curr_tracked_prns]
                df_out.at[ii, "lost_prns"] = lost_prns
                # update tracked PRNs
                tracked_prns = curr_tracked_prns
                df_out.at[ii, "num_tracked"] = len(curr_tracked_prns)

                # measurement data (take the first PRN's data as representative)
                df_out.at[ii, "cn0"] = df_tidx["cn0"].values
                df_out.at[ii, "sigma_range_m"] = df_tidx["sigma_range"].values
                df_out.at[ii, "sigma_rangerate_ms"] = df_tidx["sigma_rangerate"].values
                df_out.at[ii, "sigma_carrier_m"] = df_tidx["sigma_carrier"].values
                df_out.at[ii, "clockbias_tx_s"] = df_tidx["clockbias_tx"].values
                df_out.at[ii, "clockbias_ephem_s"] = df_tidx["clockbias_ephem"].values

                df_out.at[ii, "pos_tx_m"] = np.array(
                    [df_tidx["pos_tx"].values[i] * 1000 for i in range(len(df_tidx))]
                )
                df_out.at[ii, "vel_tx_ms"] = np.array(
                    [df_tidx["vel_tx"].values[i] * 1000 for i in range(len(df_tidx))]
                )
                df_out.at[ii, "pos_rx_m"] = df_tidx["pos_rx"].values[0] * 1000
                df_out.at[ii, "vel_rx_ms"] = df_tidx["vel_rx"].values[0] * 1000
                df_out.at[ii, "pos_tx_ephem_m"] = np.array(
                    [df_tidx["pos_tx_ephem"].values[i] * 1000 for i in range(len(df_tidx))]
                )
                df_out.at[ii, "vel_tx_ephem_ms"] = np.array(
                    [df_tidx["vel_tx_ephem"].values[i] * 1000 for i in range(len(df_tidx))]
                )

                df_out.at[ii, "total_delay_m"] = df_tidx["total_delay_m"].values
                df_out.at[ii, "tec_delay_m"] = df_tidx["tec_delay_m"].values
                df_out.at[ii, "dist_bend_m"] = df_tidx["dist_bend_m"].values
                df_out.at[ii, "second_delay_m"] = df_tidx["second_delay_m"].values
                df_out.at[ii, "third_delay_m"] = df_tidx["third_delay_m"].values
                df_out.at[ii, "total_carrier_delay_m"] = (
                    -df_tidx["tec_delay_m"].values
                    - df_tidx["tec_delay_bend_m"].values
                    - df_tidx["second_delay_m"].values / 2
                    - df_tidx["third_delay_m"].values / 3
                    + df_tidx["dist_bend_m"].values
                )
                df_out.at[ii, "tec_delay_bend_m"] = df_tidx["tec_delay_bend_m"].values
                df_out.at[ii, "min_alts_km"] = df_tidx["min_alt"].values

            # save to CSV
            output_filename = "processed_ionodata_{}_sat_{}_{}_signal_{}_rz12_{:.0f}.csv".format(
                epoch_str, lcrns_idx, gnss_const, signal_family, rz12
            )
            output_csv_path = os.path.join(savepath, output_filename)
            df_out.to_csv(output_csv_path, index=False)
            print("Saved processed data to:", output_csv_path)

            # save to pickle
            output_pickle_filename = (
                "processed_ionodata_{}_sat_{}_{}_signal_{}_rz12_{:.0f}.pkl".format(
                    epoch_str, lcrns_idx, gnss_const, signal_family, rz12
                )
            )
            output_pickle_path = os.path.join(savepath, output_pickle_filename)
            df_out.to_pickle(output_pickle_path)
            print("Saved processed data to:", output_pickle_path)

            df_outs.append(df_out)

    return df_outs


def extend_array(base_array, new_array):
    """
    Extend a dictionary by appending a new array to an existing key.
    If the key does not exist, create a new entry.
    """

    # convert to numpy array
    if base_array is not None:
        base_array = np.atleast_1d(base_array)

    if new_array is not None:
        new_array = np.atleast_1d(new_array)

    if new_array.size == 0:
        return base_array

    if base_array is not None:
        if base_array.size == 0:
            return new_array

        base_array = np.concatenate((base_array, new_array), axis=0)
    else:
        base_array = new_array

    return base_array


def add_singlefreq_measurement(
    df_tidx,
    df_meas_tidx,
    key,
    meas_rows,
    clock_bias_rx,
    tidx,
    n_mc,
    integer_ambiguity_dict,
    carrier_phase_mc_dict,
):

    gnss_const, signal_family = key.split("_")

    if signal_family == "1":
        freq_hz = 1575.42e6  # L1
        lambda_f = pnt.C / freq_hz
    elif signal_family == "5":
        freq_hz = 1176.45e6  # L5
        lambda_f = pnt.C / freq_hz
    else:
        raise ValueError(f"Unsupported signal family: {signal_family}")

    if df_tidx["tracked_prns"] is None or np.isnan(df_tidx["tracked_prns"]).all():
        return df_meas_tidx, integer_ambiguity_dict, carrier_phase_mc_dict

    # initialize measurement arrays
    prns = np.atleast_1d(df_tidx["tracked_prns"])
    acq_prns = np.atleast_1d(df_tidx["acq_prns"])
    lost_prns = np.atleast_1d(df_tidx["lost_prns"])

    n_prn = len(prns)
    t_tai = df_tidx["t_tai"]

    if len(prns) == 0:
        return df_meas_tidx, integer_ambiguity_dict, carrier_phase_mc_dict

    # extract measurement values
    total_delay = np.atleast_1d(df_tidx["total_delay_m"])  # in meters
    total_carrier_delay = np.atleast_1d(df_tidx["total_carrier_delay_m"])  # in meters
    sigma_range = np.atleast_1d(df_tidx["sigma_range_m"])  # in
    sigma_carrier = np.atleast_1d(df_tidx["sigma_carrier_m"])  # in meters
    clock_bias_tx = np.atleast_1d(df_tidx["clockbias_tx_s"])  # in seconds
    clock_bias_ephem_tx = np.atleast_1d(df_tidx["clockbias_ephem_s"])  # in seconds

    # frame conversions
    rv_tx_ecef = np.hstack([df_tidx["pos_tx_m"], df_tidx["vel_tx_ms"]])  # in meters and m/s [m, 6]
    rv_rx_ecef = np.hstack([df_tidx["pos_rx_m"], df_tidx["vel_rx_ms"]]).reshape(
        1, -1
    )  # in meters and m/s
    rv_tx_ephem_ecef = np.hstack(
        [df_tidx["pos_tx_ephem_m"], df_tidx["vel_tx_ephem_ms"]]
    )  # in meters and m/s
    rv_tx_mci = pnt.convert_frame(t_tai, rv_tx_ecef, pnt.ECEF, pnt.MOON_CI, rotate_only=False)
    rv_rx_mci = pnt.convert_frame(t_tai, rv_rx_ecef, pnt.ECEF, pnt.MOON_CI, rotate_only=False)
    rv_tx_ephem_mci = pnt.convert_frame(
        t_tai, rv_tx_ephem_ecef, pnt.ECEF, pnt.MOON_CI, rotate_only=False
    )

    range_true_m = np.linalg.norm(rv_tx_ecef[:, :3] - rv_rx_ecef[:, :3], axis=1)
    # print(f"True ranges (m): {range_true_m}")

    # construct measurements ----------------------------------------------
    for mci in range(n_mc):
        pseudorange_m = (
            range_true_m + pnt.C * (clock_bias_rx[mci, tidx, 0] - clock_bias_tx) + total_delay
        )
        carrier_phase_m = (
            range_true_m
            + pnt.C * (clock_bias_rx[mci, tidx, 0] - clock_bias_tx)
            + total_carrier_delay
        )

        # Compute integer ambiguity ----------------------------------------------
        integer_ambiguity = np.zeros(n_prn)
        for i, prn in enumerate(prns):
            if prn in df_tidx["acq_prns"]:
                # newly acquired PRN, set integer ambiguity
                integer_ambiguity[i] = np.floor(pseudorange_m[i] / lambda_f).astype(float)
                integer_ambiguity_dict[key][prn] = integer_ambiguity[i]
            elif prn in df_tidx["lost_prns"]:
                # lost PRN, remove from dictionary
                if prn in integer_ambiguity_dict[key]:
                    del integer_ambiguity_dict[key][prn]
            else:
                # existing PRN, set integer ambiguity to previous value
                integer_ambiguity[i] = integer_ambiguity_dict[key][prn]

        # print(f"Integer ambiguities for MCI {mci}, tidx {tidx}: {integer_ambiguity}")

        # add noise to measurements ----------------------------------------------
        pseudorage_noise = np.zeros(n_prn)
        carrier_phase_noise = np.zeros(n_prn)
        for i in range(n_prn):
            pseudorage_noise[i] = np.random.normal(0, sigma_range[i])
            carrier_phase_noise[i] = np.random.normal(0, sigma_carrier[i])
        pseudorange_mc = pseudorange_m + pseudorage_noise
        carrier_phase_mc = carrier_phase_m + carrier_phase_noise
        carrier_phase_phase_mc = carrier_phase_mc - integer_ambiguity * lambda_f
        graphic = (pseudorange_mc + carrier_phase_mc) / 2

        # print("Graphic values (m):", graphic)

        # register carrier phase for MC runs ----------------------------------------------
        tdcp = np.zeros(n_prn)
        for i, prn in enumerate(prns):
            if prn in df_tidx["acq_prns"]:
                # newly acquired PRN, set carrier phase
                tdcp[i] = np.nan  # no TDCP for newly acquired PRN
                carrier_phase_mc_dict[mci][key][prn] = carrier_phase_mc[i]
            elif prn in df_tidx["lost_prns"]:
                # lost PRN, remove from dictionary
                tdcp[i] = np.nan  # no TDCP for lost PRN
                if prn in carrier_phase_mc_dict[mci][key]:
                    del carrier_phase_mc_dict[mci][key][prn]
            else:
                # existing PRN, set carrier phase to previous value
                tdcp[i] = carrier_phase_mc[i] - carrier_phase_mc_dict[mci][key][prn]
                carrier_phase_mc_dict[mci][key][prn] = carrier_phase_mc[i]

        # print(f"TDCP for MCI {mci}, tidx {tidx}: {tdcp}")
        # stored measurements in dataframe ----------------------------------------------
        df_meas_tidx[mci].at[0, "tidx"] = df_tidx["tidx"]
        df_meas_tidx[mci].at[0, "tspan"] = df_tidx["tspan"]
        df_meas_tidx[mci].at[0, "t_tai"] = df_tidx["t_tai"]
        df_meas_tidx[mci].at[0, "rx_posvel_ecef"] = rv_rx_ecef
        df_meas_tidx[mci].at[0, "rx_clockstates_s"] = clock_bias_rx[mci, tidx, :]
        df_meas_tidx[mci].at[0, "rx_posvel_mci"] = rv_rx_mci

        n_acq = len(acq_prns)
        n_tracked = len(prns)
        n_lost = len(lost_prns)

        insert_vals = {
            "acq_prns": np.array(df_tidx["acq_prns"]) if n_acq > 0 else np.array([]),
            "tracked_prns": np.array(df_tidx["tracked_prns"]) if n_tracked > 0 else np.array([]),
            "lost_prns": np.array(df_tidx["lost_prns"]) if n_lost > 0 else np.array([]),
            "acq_gnss": np.array([gnss_const] * n_acq) if n_acq > 0 else np.array([]),
            "tracked_gnss": np.array([gnss_const] * n_tracked) if n_tracked > 0 else np.array([]),
            "lost_gnss": np.array([gnss_const] * n_lost) if n_lost > 0 else np.array([]),
            "acq_signals": np.array([signal_family] * n_acq) if n_acq > 0 else np.array([]),
            "tracked_signals": (
                np.array([signal_family] * n_tracked) if n_tracked > 0 else np.array([])
            ),
            "lost_signals": np.array([signal_family] * n_lost) if n_lost > 0 else np.array([]),
            "tx_posvel_ecef": rv_tx_ecef,
            "tx_clockbias_s": clock_bias_tx,
            "tx_posvel_mci": rv_tx_mci,
            "tx_posvel_ephem_ecef": rv_tx_ephem_ecef,
            "tx_clockbias_ephem_s": clock_bias_ephem_tx,
            "tx_posvel_ephem_mci": rv_tx_ephem_mci,
            "sigma_range": sigma_range,
            "sigma_carrier": sigma_carrier,
            "delay_range_m": total_delay,
            "delay_carrier_m": total_carrier_delay,
            "pseudorange": pseudorange_mc,
            "carrier_phase_m": carrier_phase_mc,
            "carrier_phase_phase": carrier_phase_phase_mc,
            "integer_ambiguity": integer_ambiguity,
            "graphic": graphic,
            "tdcp": tdcp,
        }

        for j, key_row in enumerate(meas_rows):
            if key_row not in [
                "tidx",
                "tspan",
                "t_tai",
                "rx_posvel_ecef",
                "rx_clockstates_s",
                "rx_posvel_mci",
            ] and not key_row.startswith("if_"):
                df_meas_tidx[mci].at[0, key_row] = extend_array(
                    df_meas_tidx[mci].at[0, key_row], insert_vals[key_row]
                )

    return df_meas_tidx, integer_ambiguity_dict, carrier_phase_mc_dict


def add_multifreq_measurements(df_meas_tidx, gnss_const, n_mc):
    """
    Construct Iono-free multifrequency measurements.
    """
    prns = df_meas_tidx[0].at[0, "tracked_prns"]

    if prns is None or np.isnan(prns).all():
        return df_meas_tidx

    if len(prns) == 0:
        return df_meas_tidx

    idx_gnss_mask = [i for i, key in enumerate(prns) if key == gnss_const]

    for mci in range(n_mc):
        prns = []
        ionofree_prange_list = []
        ionofree_carrier_list = []
        ionofree_sigma_range_list = []
        ionofree_sigma_carrier_list = []
        ionofree_integer_ambiguity_list = []
        ionofree_carrier_phase_phase_list = []
        if_tx_posvel_ephem_ecef_list = None
        if_tx_clockbias_ephem_s_list = None
        if_tx_posvel_ephem_mci_list = None

        for prn in np.unique(prns):
            # find indices for this PRN
            idx_prn_mask = [i for i in idx_gnss_mask if prns[i] == prn]
            if len(idx_prn_mask) < 2:
                continue  # need at least two frequencies for iono-free

            # extract measurements for the two frequencies
            prange_1 = df_meas_tidx[mci].at[0, "pseudorange"][idx_prn_mask[0]]
            prange_2 = df_meas_tidx[mci].at[0, "pseudorange"][idx_prn_mask[1]]
            carrier_1 = df_meas_tidx[mci].at[0, "carrier_phase_m"][idx_prn_mask[0]]
            carrier_2 = df_meas_tidx[mci].at[0, "carrier_phase_m"][idx_prn_mask[1]]
            sigmarange_1 = df_meas_tidx[mci].at[0, "sigma_range"][idx_prn_mask[0]]
            sigmarange_2 = df_meas_tidx[mci].at[0, "sigma_range"][idx_prn_mask[1]]
            sigmacarrier_1 = df_meas_tidx[mci].at[0, "sigma_carrier"][idx_prn_mask[0]]
            sigmacarrier_2 = df_meas_tidx[mci].at[0, "sigma_carrier"][idx_prn_mask[1]]
            iamb_1 = df_meas_tidx[mci].at[0, "integer_ambiguity"][idx_prn_mask[0]]
            iamb_2 = df_meas_tidx[mci].at[0, "integer_ambiguity"][idx_prn_mask[1]]
            cpp_1 = df_meas_tidx[mci].at[0, "carrier_phase_phase"][idx_prn_mask[0]]
            cpp_2 = df_meas_tidx[mci].at[0, "carrier_phase_phase"][idx_prn_mask[1]]
            posvel_tx_ephem_mci = df_meas_tidx[mci].at[0, "tx_posvel_ephem_mci"][idx_prn_mask[0]]
            posvel_tx_ephem_ecef = df_meas_tidx[mci].at[0, "tx_posvel_ephem_ecef"][idx_prn_mask[0]]
            clockbias_tx_ephem_s = np.array(
                [df_meas_tidx[mci].at[0, "tx_clockbias_ephem_s"][idx_prn_mask[0]]]
            )

            # frequencies
            if df_meas_tidx[mci].at[0, "tracked_signals"][idx_prn_mask[0]] == "1":
                f1 = 1575.42
                f2 = 1176.45
            else:
                f1 = 1176.45
                f2 = 1575.42

            f1_f1f2 = f1**2 / (f1**2 - f2**2)
            f2_f1f2 = f2**2 / (f1**2 - f2**2)

            ionofree_prange = f1_f1f2 * prange_1 - f2_f1f2 * prange_2
            ionofree_carrier = f1_f1f2 * carrier_1 - f2_f1f2 * carrier_2
            ionofree_sigma_range = np.sqrt(
                (f1_f1f2 * sigmarange_1) ** 2 + (f2_f1f2 * sigmarange_2) ** 2
            )
            ionofree_sigma_carrier = np.sqrt(
                (f1_f1f2 * sigmacarrier_1) ** 2 + (f2_f1f2 * sigmacarrier_2) ** 2
            )
            ionofree_integer_ambiguity = f1_f1f2 * iamb_1 - f2_f1f2 * iamb_2
            ionofree_carrier_phase_phase = f1_f1f2 * cpp_1 - f2_f1f2 * cpp_2

            ionofree_prange_list.append(ionofree_prange)
            ionofree_carrier_list.append(ionofree_carrier)
            ionofree_sigma_range_list.append(ionofree_sigma_range)
            ionofree_sigma_carrier_list.append(ionofree_sigma_carrier)
            ionofree_integer_ambiguity_list.append(ionofree_integer_ambiguity)
            ionofree_carrier_phase_phase_list.append(ionofree_carrier_phase_phase)
            if_tx_posvel_ephem_ecef_list = extend_array(
                if_tx_posvel_ephem_ecef_list, posvel_tx_ephem_ecef
            )
            if_tx_clockbias_ephem_s_list = extend_array(
                if_tx_clockbias_ephem_s_list, clockbias_tx_ephem_s
            )
            if_tx_posvel_ephem_mci_list = extend_array(
                if_tx_posvel_ephem_mci_list, posvel_tx_ephem_mci
            )
            prns.append(prn)

        # store iono-free measurements
        if len(prns) == 0:
            continue

        df_meas_tidx[mci].at[0, "if_prns"] = prns
        df_meas_tidx[mci].at[0, "if_gnss"] = [gnss_const] * len(prns)
        df_meas_tidx[mci].at[0, "if_pseudorange"] = np.array(ionofree_prange_list)
        df_meas_tidx[mci].at[0, "if_carrier_phase_m"] = np.array(ionofree_carrier_list)
        df_meas_tidx[mci].at[0, "if_sigma_range"] = np.array(ionofree_sigma_range_list)
        df_meas_tidx[mci].at[0, "if_sigma_carrier"] = np.array(ionofree_sigma_carrier_list)
        df_meas_tidx[mci].at[0, "if_integer_ambiguity"] = np.array(ionofree_integer_ambiguity_list)
        df_meas_tidx[mci].at[0, "if_carrier_phase_phase"] = np.array(
            ionofree_carrier_phase_phase_list
        )
        df_meas_tidx[mci].at[0, "if_tx_posvel_ephem_ecef"] = if_tx_posvel_ephem_ecef_list
        df_meas_tidx[mci].at[0, "if_tx_clockbias_ephem_s"] = if_tx_clockbias_ephem_s_list
        df_meas_tidx[mci].at[0, "if_tx_posvel_ephem_mci"] = if_tx_posvel_ephem_mci_list

    return df_meas_tidx


def generate_iono_gnss_meas(sim_params, clock_model="OCXO", n_mc=10, seed=0):
    """
    Generate GNSS measurement data from processed ionodata CSV files.
    """
    epoch_str = sim_params["epoch_str"]
    lcrns_idx = sim_params["lcrns_idx"]
    rz12 = sim_params["rz12"]
    n_orbit = sim_params["n_orbit"]
    dt = sim_params["dt"]
    dt_raytrace = sim_params["dt_raytrace"]
    df_dict = {}

    orbit_dt_str = "norbit_{}_dt_{}s_dtrt_{}s".format(n_orbit, dt, dt_raytrace)

    datapath = os.path.join(pnt.get_output_dir(), "iono_delay", "raytrace_timestep", orbit_dt_str)

    t0 = np.inf
    tf = -np.inf
    dt = np.inf

    # load the pre-generated pickle files
    for gnss_const in sim_params["gnss_consts"]:
        for signal_family in sim_params["signal_families"]:
            filepath = os.path.join(
                datapath,
                "processed_ionodata_{}_sat_{}_{}_signal_{}_rz12_{:.0f}.pkl".format(
                    epoch_str, lcrns_idx, gnss_const, signal_family, rz12
                ),
            )
            df = pd.read_pickle(filepath)
            print("Loaded GNSS measurement data from:", filepath)
            df_dict[f"{gnss_const}_{signal_family}"] = df

            tspan = df["tspan"].values
            t0 = min(tspan[0], t0)
            tf = max(tspan[-1], tf)
            dt = min(tspan[1] - tspan[0], dt)

    # construct global time span
    tspan = np.arange(t0, tf + dt, dt)
    lent = len(tspan)

    # simulate clocks --------------------------------------------------------------
    print("Simulating receiver clock bias...")
    clock_noise = ClockNoise(clock_model)
    clock_bias_rx = np.zeros((n_mc, len(tspan), 3))
    for mci in range(n_mc):
        clock_bias_rx[mci, :, :] = clock_noise.simulate_clock_bias(tspan)

    # construct measurements
    meas_rows = [
        "tidx",
        "tspan",
        "t_tai",
        "acq_prns",
        "tracked_prns",
        "lost_prns",
        "acq_gnss",
        "tracked_gnss",
        "lost_gnss",
        "acq_signals",
        "tracked_signals",
        "lost_signals",
        "rx_posvel_ecef",
        "rx_clockstates_s",
        "rx_posvel_mci",
        "tx_posvel_ecef",
        "tx_clockbias_s",
        "tx_posvel_mci",
        "tx_posvel_ephem_ecef",
        "tx_clockbias_ephem_s",
        "tx_posvel_ephem_mci",
        "sigma_range",
        "sigma_carrier",
        "delay_range_m",
        "delay_carrier_m",
        "pseudorange",
        "carrier_phase_m",
        "carrier_phase_phase",
        "integer_ambiguity",
        "graphic",
        "tdcp",
        # Iono-free measurements
        "if_prns",
        "if_gnss",
        "if_pseudorange",
        "if_carrier_phase_m",
        "if_sigma_range",
        "if_sigma_carrier",
        "if_integer_ambiguity",
        "if_carrier_phase_phase",
        "if_tx_posvel_ephem_ecef",
        "if_tx_clockbias_ephem_s",
        "if_tx_posvel_ephem_mci",
    ]
    df_meas = pd.DataFrame(columns=meas_rows, index=np.arange(lent))
    df_meas_mc = [df_meas.copy() for _ in range(n_mc)]

    # Generate Measurements --------------------------------------------------------------
    # Todo:
    # - GPS-GALILEO inter-frequency biases
    # - relativistic corrections

    # Initialize
    integer_ambiguity_dict = {}
    carrier_phase_mc_dict = [{} for _ in range(n_mc)]
    for key in df_dict.keys():
        integer_ambiguity_dict[key] = {}
        for mci in range(n_mc):
            carrier_phase_mc_dict[mci][key] = {}

    np.random.seed(seed)

    for tidx in tqdm(range(lent), desc="Generating GNSS Measurements"):
        # Initialize measurement rows for this time index ----------------------------------
        for meas_col in meas_rows:
            for mci in range(n_mc):
                df_meas_mc[mci].at[tidx, meas_col] = None

        df_meas_tidx = [pd.DataFrame(columns=meas_rows, index=[0]) for _ in range(n_mc)]
        for meas_col in meas_rows:
            for mci in range(n_mc):
                df_meas_tidx[mci].at[0, meas_col] = None

        # First construct measurement for each gnss and frequencies ----------------------------------
        # print("  Adding single-frequency measurements...")
        for key in df_dict.keys():
            df_tidx = df_dict[key].iloc[tidx, :]
            gnss_const, signal_family = key.split("_")

            df_meas_tidx, integer_ambiguity_dict, carrier_phase_mc_dict = (
                add_singlefreq_measurement(
                    df_tidx,
                    df_meas_tidx,
                    key,
                    meas_rows,
                    clock_bias_rx,
                    tidx,
                    n_mc,
                    integer_ambiguity_dict,
                    carrier_phase_mc_dict,
                )
            )

        # Now aggregate all measurements across frequencies ----------------------------------===========
        # print("  Adding multifrequency measurements...")
        for gnss_const in sim_params["gnss_consts"]:
            # concatenate measurements across frequencies
            df_meas_tidx = add_multifreq_measurements(df_meas_tidx, gnss_const, n_mc)

        # Now store the measurements for all MC runs ----------------------------------===========
        # print("  Storing measurements for all MC runs...")
        for mci in range(n_mc):
            for col in meas_rows:
                df_meas_mc[mci].at[tidx, col] = df_meas_tidx[mci].at[0, col]

        # print("Completed storing measurements for time index:", tidx)

    # save the df for each mc run --------------------------------------------------------------
    savepath = os.path.join(pnt.get_output_dir(), "iono_delay", "gnss_measurements", orbit_dt_str)
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    for mci in range(n_mc):
        output_filename = "gnss_measurements_iono_{}_sat_{}_rz12_{:.0f}_mc_{:02d}.pkl".format(
            epoch_str, lcrns_idx, rz12, mci
        )
        output_filepath = os.path.join(savepath, output_filename)
        df_meas_mc[mci].to_pickle(output_filepath)
        print("Saved GNSS measurement data to:", output_filepath)

    return df_meas_mc
