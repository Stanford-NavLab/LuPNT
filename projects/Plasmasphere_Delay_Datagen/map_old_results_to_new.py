"""
Function to map old raytrace results from the old raytrace CSV files into the new dataframe. It handles both normal
and corrected raytrace results
"""

import signal
from unittest import signals
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import pylupnt as pnt


def load_files(sim_params, datadir_labels_new, datadir_labels_old, datadir_rt_results_old):
    ymdh = sim_params["epoch_ymdh"]
    lcrns_idx = sim_params["lcrns_idx"]
    rz12 = sim_params["rz12"]
    gnss_const = sim_params["gnss_const"]
    signal = sim_params["signal"]
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

    # load new labels Pickle ---------------------------------------------------------------
    pickle_labels_filename_new = os.path.join(
        datadir_labels_new,
        "ionodata_short_{}_sat_{}_{}_signal_{}.pkl".format(
            epoch_str, lcrns_idx, gnss_const, signal
        ),
    )

    if not os.path.exists(pickle_labels_filename_new):
        raise FileNotFoundError(f"Pickle file {pickle_labels_filename_new} does not exist.")

    # Load the initial ionodata Pickle
    df_new = pd.read_pickle(pickle_labels_filename_new)
    print(f"Pickle file {pickle_labels_filename_new} loaded successfully.")
    print("Number of rows in the new dataframe:", len(df_new))

    # Load the labels CSV for the old one ---------------------------------------------------------------
    pickle_labels_filename_old = os.path.join(
        datadir_labels_old,
        "ionodata_short_{}_sat_{}_{}_signal_{}.pkl".format(
            epoch_str, lcrns_idx, gnss_const, signal
        ),
    )
    if not os.path.exists(pickle_labels_filename_old):
        raise FileNotFoundError(f"Pickle file {pickle_labels_filename_old} does not exist.")

    df_old = pd.read_pickle(pickle_labels_filename_old)
    print(f"Pickle file {pickle_labels_filename_old} loaded successfully.")
    print("Number of rows in the old dataframe:", len(df_old))

    # Load the delay data CSV for the old one ---------------------------------------------------------------
    ray_trace_csv_dir_old = os.path.join(
        datadir_rt_results_old,
        "raytrace_results_sat_{}_{}_signal_{}_rz12_{:.0f}_kp_{:.0f}".format(
            lcrns_idx, gnss_const, signal, rz12, kp
        ),
    )

    ray_trace_csv_dir_old_correction = os.path.join(
        datadir_rt_results_old,
        "raytrace_results_correct_sat_{}_{}_signal_{}_rz12_{:.0f}_kp_{:.0f}".format(
            lcrns_idx, gnss_const, signal, rz12, kp
        ),
    )

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

    # Load and overwrite delay data from ray trace CSV files (normal)
    ray_trace_files = [f for f in os.listdir(ray_trace_csv_dir_old) if f.endswith(".csv")]
    print(
        "Number of ray trace CSV files found in ", ray_trace_csv_dir_old, ":", len(ray_trace_files)
    )
    for file in tqdm(ray_trace_files, desc="Loading ray trace CSV files"):
        # extract row number from filename (raytrace_result_{row}.csv)
        row_num = int(file.split("_")[-1].split(".")[0])
        df_ray = pd.read_csv(os.path.join(ray_trace_csv_dir_old, file))
        # replace columns
        # total_delay_m,tecu,tec_delay_m,second_delay_m,third_delay_m,dist_bend_m,tec_delay_bend_m,max_sep_line_m,final_pos_err_m
        data = df_ray[raytrace_cols].iloc[0]
        # replace the corresponding row in df
        for col in raytrace_cols:
            df_old.at[row_num, col] = data[col]

    # Load and overwrite delay data from ray trace CSV files (correction)
    ray_trace_files_correction = [
        f for f in os.listdir(ray_trace_csv_dir_old_correction) if f.endswith(".csv")
    ]
    print(
        "Number of ray trace CSV files found in ",
        ray_trace_csv_dir_old_correction,
        ":",
        len(ray_trace_files_correction),
    )
    for file in tqdm(ray_trace_files_correction, desc="Loading ray trace CSV files (correction)"):
        # extract row number from filename (raytrace_result_{row}.csv)
        row_num = int(file.split("_")[-1].split(".")[0])
        df_ray_correct = pd.read_csv(os.path.join(ray_trace_csv_dir_old_correction, file))
        # replace columns
        data = df_ray_correct[raytrace_cols].iloc[0]
        # replace the corresponding row in df
        for col in raytrace_cols:
            df_old.at[row_num, col] = data[col]

    return df_old, df_new


def convert_old_to_new(sim_params, df_old, df_new, datadir_rt_results_new, overwrite=False):
    """Convert old ephemeris data to new ephemeris data using scaling factors.
    Args:
        sim_params (dict): Simulation parameters containing 'gnss_const', 'lcrns_idx', 'epoch_ymdh', and 'rz12'.
        df_old (pd.DataFrame): DataFrame containing the old ephemeris data.
        df_new (pd.DataFrame): DataFrame containing the new ephemeris data.
        datadir_rt_results_new (str): Directory path where the new ray trace CSV files are stored.
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

    datadir_save = os.path.join(
        datadir_rt_results_new,
        "raytrace_results_sat_{}_{}_signal_{}_rz12_{:.0f}_kp_{:.0f}".format(
            sim_params["lcrns_idx"],
            sim_params["gnss_const"],
            sim_params["signal"],
            sim_params["rz12"],
            sim_params["kp"],
        ),
    )
    datadir_correct_save = os.path.join(
        datadir_rt_results_new,
        "raytrace_results_correct_sat_{}_{}_signal_{}_rz12_{:.0f}_kp_{:.0f}".format(
            sim_params["lcrns_idx"],
            sim_params["gnss_const"],
            sim_params["signal"],
            sim_params["rz12"],
            sim_params["kp"],
        ),
    )

    os.makedirs(datadir_save, exist_ok=True)
    os.makedirs(datadir_correct_save, exist_ok=True)

    # Iterate over each row in df_new and find the corresponding L1/E1 data
    L1_tidx_row = 0

    for idx, row in tqdm(
        df_new.iterrows(), total=len(df_new), desc="Creating new data with old data"
    ):
        tidx = row["tidx"]

        # Update L1_tidx_row_start and L1_tidx_row_end to narrow down the search rang
        old_rows_tidx = []
        for old_tidx_row in range(L1_tidx_row, len(df_old)):
            if df_old.at[old_tidx_row, "tidx"] == tidx:
                old_rows_tidx.append(old_tidx_row)
            if old_tidx_row >= len(df_old):
                break
            if df_old.at[old_tidx_row, "tidx"] > tidx:
                break

        # Find the matching row in df_old based on min_alt
        match_old = df_old.loc[old_rows_tidx]
        match_old = match_old[(match_old["min_alt"] == row["min_alt"])]

        if not match_old.empty:
            # There should be only one matching row
            df_data = match_old.iloc[0]
            # Copy ephemeris data from match_row to the current row in df
            df_L5_row = pd.DataFrame(
                {
                    "total_delay_m": df_data["total_delay_m"],
                    "tecu": df_data["tecu"],
                    "tec_delay_m": df_data["tec_delay_m"],
                    "second_delay_m": df_data["second_delay_m"],
                    "third_delay_m": df_data["third_delay_m"],
                    "dist_bend_m": df_data["dist_bend_m"],
                    "tec_delay_bend_m": df_data["tec_delay_bend_m"],
                    "max_sep_line_m": df_data["max_sep_line_m"],
                    "final_pos_err_m": df_data["final_pos_err_m"],
                },
                index=[0],
            )
            match_num += 1

            # save the updated df_new to CSV
            df_rt_new_filename = os.path.join(datadir_save, "raytrace_result_{}.csv".format(idx))

            # save the updated df_new to CSV
            if not os.path.exists(df_rt_new_filename) or overwrite:
                # save only if file does not exist
                df_L5_row.to_csv(df_rt_new_filename, index=False)

            if df_data["dist_bend_m"] > 0:
                # raytrace with bending corrections
                df_rt_correct_filename = os.path.join(
                    datadir_correct_save, "raytrace_result_{}.csv".format(idx)
                )
                if not os.path.exists(df_rt_correct_filename) or overwrite:
                    df_L5_row.to_csv(df_rt_correct_filename, index=False)

    # print how many L5/E5a rows were matched with L1/E1 data
    print(f"Number of matched L5/E5a rows with L1/E1 data: {match_num}/{len(df_new)}")


if __name__ == "__main__":

    gnss_consts = ["GPS"]
    signal_families = [1, 5]
    lcrns_idx = 0
    epoch_ymdh = [2025, 3, 1, 12]
    rz12 = 50.0
    kp = 3.0

    dt = 1
    dtrt = 120
    n_orbit = 6

    datapath_labels_old = os.path.join(
        pnt.get_output_dir(),
        "iono_delay",
        "labels_csv",
        f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s_gpsold",
    )
    datapath_raytrace_old = os.path.join(
        pnt.get_output_dir(),
        "iono_delay",
        "raytrace_csv",
        f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s_gpsold",
    )
    datapath_labels_new = os.path.join(
        pnt.get_output_dir(), "iono_delay", "labels_csv", f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s"
    )
    datapath_raytrace_new = os.path.join(
        pnt.get_output_dir(),
        "iono_delay",
        "raytrace_csv",
        f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s",
    )

    for gnss_const in gnss_consts:
        for signal_family in signal_families:

            print("--------------------------------------------")
            print(f"Processing GNSS constellation: {gnss_const}, Signal family: {signal_family}")
            print(f"Datapath Old Labels: {datapath_labels_old}")
            print(f"Datapath Old Raytrace: {datapath_raytrace_old}")
            print(f"Datapath New Labels: {datapath_labels_new}")
            print(f"Datapath New Raytrace: {datapath_raytrace_new}")
            print("--------------------------------------------")

            sim_params = {
                "gnss_const": gnss_const,
                "lcrns_idx": lcrns_idx,
                "epoch_ymdh": epoch_ymdh,
                "rz12": rz12,  # -1.0 for historical or projected R12
                "kp": kp,
                "signal": signal_family,
            }

            df_old, df_new = load_files(
                sim_params, datapath_labels_new, datapath_labels_old, datapath_raytrace_old
            )
            convert_old_to_new(sim_params, df_old, df_new, datapath_raytrace_new, overwrite=False)
