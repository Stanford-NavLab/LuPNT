import numpy as np
import pandas as pd
import pylupnt as pnt
from tqdm import tqdm
from datetime import datetime, timedelta
from src.setup_ionoenv import setup_lcrns_sats, setup_gnss_constellation
from src.generate_ionodata import generate_labels

# setting
gps_datetime = datetime(2025, 3, 1, 12, 0, 0)  # set the GPS datetime for the simulation
n_orbit = 6.0  # number of orbits to simulate for TDCP
dt = 1.0  # time step in seconds
tidx_inv_raytrace = 120  # time index for inverse raytracing (6 corresponds to 1 minute)
overwrite_rx_orbit = False  # whether to overwrite existing orbit data
overwrite_gnss_orbit = False  # whether to overwrite existing GNSS data
overwrite_gnss_measurements = True  # whether to overwrite existing GNSS measurements
overwrite_labels = True  # whether to overwrite existing ionosphere label data
consider_gnss_fault = True  # whether to consider GNSS faults
sp3_prop_method = "interp"  # method for SP3 propagation ("interp" gives better accuracy)

# overwrite coupled paraemters
if overwrite_rx_orbit:
    overwrite_gnss_measurements = True
    overwrite_labels = True
elif overwrite_gnss_orbit:
    overwrite_gnss_measurements = True
    overwrite_labels = True
elif overwrite_gnss_measurements:
    overwrite_labels = True

# Setup the simulation settings
t_tai, rv_m2sc_ci, rv_e2sc_ecef, N_sc = setup_lcrns_sats(
    n_orbit=n_orbit, dt=dt, savefig=False, overwrite=overwrite_rx_orbit
)

# setup GNSS constellation and measurements
gnss_meas = setup_gnss_constellation(
    t_tai,
    rv_m2sc_ci,
    rv_e2sc_ecef=None,
    gps_datetime=gps_datetime,
    savefig=False,
    overwrite_orbit=overwrite_gnss_orbit,
    overwrite_measurements=overwrite_gnss_measurements,
    sp3_prop_method=sp3_prop_method,
    consider_faults=consider_gnss_fault,
)
N_t = gnss_meas.N_t  # number of time steps

# Generate labels for different configurations ------------------------------------------------------------------
signal_families = [1, 5]
gnss_consts = ["GPS"]  # , "GALILEO", "QZSS"]
lcrns_idxs = [0]  # LCRNS 1 and South Pole

for signal_family in signal_families:
    for gnss_const in gnss_consts:
        for lcrns_idx in lcrns_idxs:
            print(
                f"Generating labels for {gnss_const} signal family {signal_family} with LCRNS index {lcrns_idx}..."
            )
            sim_params = {
                "gnss_const": gnss_const,
                "signal_family": signal_family,
                "lcrns_idx": lcrns_idx,
                "epoch_ymdh": [2025, 3, 1, 12],
                "correction": False,
                "dt": dt,
                "use_lt": True,
                "n_orbit": n_orbit,
            }
            generate_labels(
                gnss_meas,
                sim_params,
                overwrite=overwrite_labels,
                tidxs_inv_raytrace=tidx_inv_raytrace,
            )
