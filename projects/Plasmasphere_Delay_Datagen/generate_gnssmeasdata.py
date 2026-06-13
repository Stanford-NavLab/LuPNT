import pylupnt as pnt
import pandas as pd
import numpy as np
import os
from src.postprocess_ionodata import generate_iono_gnss_meas

# simulation parameters
sim_params = {
    "gnss_consts": ["GPS", "GALILEO"],
    "signal_families": [1, 5],
    "lcrns_idx": 0,
    "epoch_str": "2025_03_01_12_00_00",
    "rz12": 50.0,
    "kp": 3.0,
    "n_orbit": 3,
    "dt": 1,  # time step in seconds
    "dt_raytrace": 60,  # raytracing time step in seconds
}

n_mc = 10

meas_dict_list = generate_iono_gnss_meas(
    sim_params=sim_params, clock_model="OCXO", n_mc=n_mc, seed=0
)
