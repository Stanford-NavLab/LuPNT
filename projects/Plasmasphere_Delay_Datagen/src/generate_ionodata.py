import numpy as np
import pandas as pd
from tqdm import tqdm
import tecsimpy as tec
import pylupnt as pnt
import os


def generate_labels(gnss_meas, sim_params, overwrite=False, tidxs_inv_raytrace=6):
    """
    Generate ionosphere delay labels for the given GNSS measurements and simulation parameters.

    Parameters
    ----------
    gnss_meas : pnt.GNSSMeas
        GNSS measurement object containing visibility and measurement data.
    sim_params : dict
        Simulation parameters including:
    overwrite : bool, optional
        Whether to overwrite existing label files, by default False.
    tidxs_inv_raytrace : int, optional
        Interval for raytrace labels, by default 6. (every 60 seconds if dt=10s)
    """

    gnss_const = sim_params["gnss_const"]
    signal_family = sim_params["signal_family"]
    lcrns_idx = sim_params["lcrns_idx"]
    min_cn0 = sim_params.get("min_cn0", 18.0)
    dt = sim_params["dt"]
    use_lt = sim_params.get("use_lt", True)
    n_orbit = sim_params.get("n_orbit", 3.0)

    N_t = gnss_meas.N_t  # number of time steps
    tspan = gnss_meas.tspan  # time span in seconds since first epoch

    # time settings ------------------------------------------------
    epoch_dict = {
        "year": sim_params["epoch_ymdh"][0],
        "month": sim_params["epoch_ymdh"][1],
        "day": sim_params["epoch_ymdh"][2],
        "hour": sim_params["epoch_ymdh"][3],
        "minute": 0,
        "second": 0,
    }
    epoch_str = "{year}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}_{second:02d}".format(
        **epoch_dict
    )
    mjd = pnt.gregorian_to_mjd(
        epoch_dict["year"],
        epoch_dict["month"],
        epoch_dict["day"],
        epoch_dict["hour"],
        epoch_dict["minute"],
        epoch_dict["second"],
    )
    epoch_utc = tec.mjd_to_tj2000(mjd)

    print("dt:", dt)
    tidxs = np.arange(N_t)
    rt_start_idx = 1  # first measurement starts from t=1

    tidxs_rt = np.arange(rt_start_idx, N_t, tidxs_inv_raytrace)
    print("tidx:", tidxs)
    print("tidx raytrace:", tidxs_rt)

    # create signal name -------------------------------------------
    if signal_family == 1:
        signal = "L1" if (gnss_const == "GPS" or gnss_const == "QZSS") else "E1"
    elif signal_family == 2:
        signal = "L2" if (gnss_const == "GPS" or gnss_const == "QZSS") else "E6"
    elif signal_family == 5:
        signal = "L5" if (gnss_const == "GPS" or gnss_const == "QZSS") else "E5a"
    else:
        raise ValueError(f"Unsupported signal family: {signal_family}")

    # check if labels file exists --------------------------------
    filename = f"ionodata_full_{epoch_str}_sat_{lcrns_idx}_{gnss_const}_signal_{signal_family}.pkl"
    orbit_str = "norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s".format(
        n_orbit=int(n_orbit), dt=int(dt), dtrt=int(dt * tidxs_inv_raytrace)
    )
    filedir = os.path.join(pnt.get_output_dir(), "iono_delay", "labels_csv", orbit_str)
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    filename_labels = os.path.join(filedir, filename)

    filename_rt_pkl = (
        f"ionodata_short_{epoch_str}_sat_{lcrns_idx}_{gnss_const}_signal_{signal_family}.pkl"
    )
    filename_labels_rt = os.path.join(filedir, filename_rt_pkl)

    df_exists = False
    df_rt_exists = False
    if filename_labels is not None and os.path.exists(filename_labels) and not overwrite:
        df_exists = True
    if filename_labels_rt is not None and os.path.exists(filename_labels_rt) and not overwrite:
        df_rt_exists = True

    if df_exists and df_rt_exists:
        print("    Labels already exist, skipping generation...")
        return

    print("    Generating new labels...")

    # Create CSV with obtained data ------------------------------
    table_rows = [
        "tidx",
        "t_tai",
        "tspan",
        "epoch_t",
        "gnss_const",
        "signal",
        "sat_id",
        "prn",
        "min_alt",
        "sigma_range",
        "sigma_rangerate",
        "sigma_carrier",
        "cn0",
        "pos_tx",
        "pos_rx",
        "vel_tx",
        "vel_rx",
        "clockbias_tx",
        "pos_tx_ephem",
        "vel_tx_ephem",
        "clockbias_ephem",
        "tecu",
        "tec_delay_m",
        "second_delay_m",
        "third_delay_m",
        "dist_bend_m",
        "tec_delay_bend_m",
        "total_delay_m",
        "max_sep_line_m",
        "final_pos_err_m",
    ]

    table_rows_rt = [
        "tidx",
        "row_full",
        "epoch_t",
        "min_alt",
        "pos_tx",
        "pos_rx",
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
    # first count the number of valid labels
    num_labels = 0
    num_labels_rt = 0
    for tidx in tidxs:
        t = tspan[tidx]  # time since first epoch in seconds
        vis_t = gnss_meas.vis_gnss[gnss_const][signal][lcrns_idx, :, tidx]
        gps_idxs = np.where(vis_t)[0]

        for gps_idx in gps_idxs:
            cn0 = gnss_meas.cn0[gnss_const][signal][lcrns_idx, gps_idx, tidx]
            if cn0 < min_cn0:  # skip low C/N0 measurements
                continue
            num_labels += 1
            if tidx in tidxs_rt:
                num_labels_rt += 1

    df = pd.DataFrame(columns=table_rows, index=np.arange(num_labels))
    df_raytrace = pd.DataFrame(columns=table_rows_rt, index=np.arange(num_labels_rt))
    print(f"    Number of labels to be generated: {num_labels}")
    print(f"    Number of raytrace labels to be generated: {num_labels_rt}")

    label_dict_idx = 0
    short_idx = 0

    # if not existing, generate new labels
    for tidx in tqdm(tidxs, desc="Generating labels"):
        t = tspan[tidx]  # time since first epoch in seconds
        vis_t = gnss_meas.vis_gnss[gnss_const][signal][lcrns_idx, :, tidx]
        gps_idxs = np.where(vis_t)[0]

        # extract ionospheric parameters
        # ionoparams = tec.get_iono_params(epoch_utc + t, kp)

        for gps_idx in gps_idxs:
            cn0 = gnss_meas.cn0[gnss_const][signal][lcrns_idx, gps_idx, tidx]
            if cn0 < min_cn0:  # skip low C/N0 measurements
                continue

            pos_tx = gnss_meas.rv_gnss_ecef[gnss_const][
                gps_idx, tidx, :3
            ]  # position of the GNSS satellite in CI frame [km]
            pos_rx = gnss_meas.rv_e2sc_ecef[
                lcrns_idx, tidx, :3
            ]  # position of the LCRNS in CI frame [km]
            vel_rx = gnss_meas.rv_e2sc_ecef[
                lcrns_idx, tidx, 3:
            ]  # velocity of the LCRNS in CI frame [km/s]

            # compute epoch
            if use_lt:
                rv_tx_new, lt = gnss_meas.lt_correction(gnss_const, lcrns_idx, gps_idx, tidx)
                pos_tx_new = rv_tx_new[:3]
                t_tx = t - lt  # apply light time correction
                pos_tx = pos_tx_new  # update position of the GNSS satellite
            else:
                rv_tx_new = gnss_meas.rv_gnss_ecef[gnss_const][gps_idx, tidx, :]
                t_tx = t

            # covert to km, km/s
            pos_tx = rv_tx_new[:3] * 1e-3  # convert to km
            vel_tx = rv_tx_new[3:] * 1e-3  # convert to km/s
            pos_rx = pos_rx * 1e-3  # convert to km
            vel_rx = vel_rx * 1e-3  # convert to km/s
            sigma_range = (
                gnss_meas.sigma_range[gnss_const][signal][lcrns_idx, gps_idx, tidx] * 1e-3
            )  # convert to km
            sigma_rangerate = (
                gnss_meas.sigma_rangerate[gnss_const][signal][lcrns_idx, gps_idx, tidx] * 1e-3
            )  # convert to km/s
            sigma_carrier_phase = (
                gnss_meas.sigma_carrier_phase[gnss_const][signal][lcrns_idx, gps_idx, tidx] * 1e-3
            )  # convert to km

            epoch_t = epoch_utc + t_tx  # epoch for the current time step

            min_alt = tec.compute_min_altitude(pos_tx, pos_rx, tec.RE)

            if min_alt <= 100.0:
                # skip if min altitude is below 100 km
                continue

            if gnss_const == "GPS":
                gnss_str = "G"
                prn = gnss_meas.prn_gnss[gnss_const][gps_idx]
            elif gnss_const == "GALILEO":
                gnss_str = "E"
                prn = gnss_meas.prn_gnss[gnss_const][gps_idx]
            elif gnss_const == "QZSS":
                gnss_str = "J"
                prn = gnss_meas.prn_gnss[gnss_const][gps_idx]
            else:
                raise ValueError(f"Unsupported GNSS constellation: {gnss_const}")

            t_tai_tx = gnss_meas.t_tai[tidx] - lt
            rv_true, clock_true = gnss_meas.sp3l.get_posvelclock(
                gnss_str, prn, t_tai_tx, out_frame=pnt.ECEF
            )
            rv_ephem, clock_ephem = gnss_meas.brdc.get_posvelclock(
                gnss_str, prn, t_tai_tx, out_frame=pnt.ECEF
            )
            rv_diff = rv_ephem - rv_true
            pos_tx_ephem = pos_tx + rv_diff[:3] * 1e-3  # convert to km
            vel_tx_ephem = vel_tx + rv_diff[3:] * 1e-3  # convert to km/s

            # full dataframe
            df.at[label_dict_idx, "tidx"] = tidx
            df.at[label_dict_idx, "t_tai"] = gnss_meas.t_tai[tidx]
            df.at[label_dict_idx, "tspan"] = gnss_meas.tspan[tidx]
            df.at[label_dict_idx, "epoch_t"] = epoch_t
            df.at[label_dict_idx, "gnss_const"] = gnss_const
            df.at[label_dict_idx, "signal"] = signal
            df.at[label_dict_idx, "sat_id"] = lcrns_idx
            df.at[label_dict_idx, "prn"] = prn
            df.at[label_dict_idx, "min_alt"] = min_alt
            df.at[label_dict_idx, "sigma_range"] = sigma_range
            df.at[label_dict_idx, "sigma_rangerate"] = sigma_rangerate
            df.at[label_dict_idx, "sigma_carrier"] = sigma_carrier_phase
            df.at[label_dict_idx, "cn0"] = cn0
            df.at[label_dict_idx, "pos_tx"] = pos_tx
            df.at[label_dict_idx, "pos_rx"] = pos_rx
            df.at[label_dict_idx, "vel_tx"] = vel_tx
            df.at[label_dict_idx, "vel_rx"] = vel_rx
            df.at[label_dict_idx, "clockbias_tx"] = clock_true
            df.at[label_dict_idx, "pos_tx_ephem"] = pos_tx_ephem
            df.at[label_dict_idx, "vel_tx_ephem"] = vel_tx_ephem
            df.at[label_dict_idx, "clockbias_ephem"] = clock_ephem
            # set the remaining labels to zero for now
            df.at[label_dict_idx, "tecu"] = 0.0
            df.at[label_dict_idx, "tec_delay_m"] = 0.0
            df.at[label_dict_idx, "second_delay_m"] = 0.0
            df.at[label_dict_idx, "third_delay_m"] = 0.0
            df.at[label_dict_idx, "dist_bend_m"] = 0.0
            df.at[label_dict_idx, "tec_delay_bend_m"] = 0.0
            df.at[label_dict_idx, "total_delay_m"] = 0.0
            df.at[label_dict_idx, "max_sep_line_m"] = 0.0
            df.at[label_dict_idx, "final_pos_err_m"] = 0.0

            # add to the short list if tidx is in tidxs_short
            if tidx in tidxs_rt:
                df_raytrace.at[short_idx, "tidx"] = df.at[label_dict_idx, "tidx"]
                df_raytrace.at[short_idx, "row_full"] = label_dict_idx
                df_raytrace.at[short_idx, "epoch_t"] = df.at[label_dict_idx, "epoch_t"]
                df_raytrace.at[short_idx, "min_alt"] = df.at[label_dict_idx, "min_alt"]
                df_raytrace.at[short_idx, "pos_tx"] = df.at[label_dict_idx, "pos_tx"]
                df_raytrace.at[short_idx, "pos_rx"] = df.at[label_dict_idx, "pos_rx"]
                df_raytrace.at[short_idx, "total_delay_m"] = 0.0
                df_raytrace.at[short_idx, "tecu"] = 0.0
                df_raytrace.at[short_idx, "tec_delay_m"] = 0.0
                df_raytrace.at[short_idx, "second_delay_m"] = 0.0
                df_raytrace.at[short_idx, "third_delay_m"] = 0.0
                df_raytrace.at[short_idx, "dist_bend_m"] = 0.0
                df_raytrace.at[short_idx, "tec_delay_bend_m"] = 0.0
                df_raytrace.at[short_idx, "max_sep_line_m"] = 0.0
                df_raytrace.at[short_idx, "final_pos_err_m"] = 0.0
                short_idx += 1

            # increment counter
            label_dict_idx += 1

    # extract only the filled rows
    df = df[:label_dict_idx]
    df_raytrace = df_raytrace[:short_idx]

    num_labels = len(df)
    print(f"    Number of labels generated: {num_labels}")
    num_raytrace_labels = len(df_raytrace)
    print(f"    Number of raytrace labels generated: {num_raytrace_labels}")

    # save to file ---------------------------------------------
    if not df_exists:
        df.to_pickle(filename_labels)
        # df.to_csv(filename_labels_csv, index=False)
    if not df_rt_exists:
        df_raytrace.to_pickle(filename_labels_rt)
        # df_raytrace.to_csv(filename_labels_rt_csv, index=False)

    # return df, df_raytrace
