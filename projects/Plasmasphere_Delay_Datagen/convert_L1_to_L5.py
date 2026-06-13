import os
import numpy as np
import pandas as pd
import pylupnt as pnt
from tqdm import tqdm
from src.postprocess_ionodata import load_L1_L5_labels, convert_L1_to_L5

# Parameters --------------------------
epoch_ymdh = [2025, 3, 1, 12]
lcrns_idxs = [0]
rz12s = [50.0]
gnss_consts = ["GPS", "GALILEO", "QZSS"]
kp = 3.0

n_orbit = 6
dt = 1
dtrt = 120
datapath_labels = os.path.join(
    pnt.get_output_dir(), "iono_delay", "labels_csv", f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s"
)
datapath_raytrace = os.path.join(
    pnt.get_output_dir(), "iono_delay", "raytrace_csv", f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s"
)
overwrite = True

# conversions
for lcrn_idx in lcrns_idxs:
    for rz12 in rz12s:
        for gnss_const in gnss_consts:
            sim_params = {
                "gnss_const": gnss_const,
                "lcrns_idx": lcrn_idx,
                "epoch_ymdh": epoch_ymdh,
                "rz12": rz12,  # -1.0 for historical or projected R12
                "kp": kp,
            }

            print(f"Processing GNSS const: {gnss_const}, LCRNS idx: {lcrn_idx}, RZ12: {rz12}...")

            df_L1, df_L5 = load_L1_L5_labels(sim_params, datapath_labels, datapath_raytrace)
            df_L5_converted = convert_L1_to_L5(
                sim_params, df_L1, df_L5, datapath_raytrace, overwrite
            )

            print(
                f"Converted L1 to L5 for GNSS const: {gnss_const}, LCRNS idx: {lcrn_idx}, RZ12: {rz12}"
            )
