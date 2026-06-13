import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import pylupnt as pnt
from src.postprocess_ionodata import generate_meas_timestep_df

# datapaths
datapath = os.path.join(pnt.get_output_dir(), "iono_delay", "raytrace_summary")
savepath = os.path.join(pnt.get_output_dir(), "iono_delay", "raytrace_timestep")

# list of all combinations to run
COMBOS = [
    ("GPS", 1),
    ("GPS", 5),
    # ("GALILEO", 1),
    # ("GALILEO", 5),
]


def run_single(combo, worker_id=0):
    """Worker function for one (gnss_const, signal_family) pair."""
    gnss_const, signal_family = combo

    sim_params = {
        "gnss_consts": [gnss_const],
        "signal_families": [signal_family],
        "lcrns_idx": 0,
        "epoch_str": "2025_03_01_12_00_00",
        "rz12": 50.0,
        "kp": 3.0,
        "n_orbit": 3,
        "dt": 1,
        "dt_raytrace": 60,
    }

    df_outs = generate_meas_timestep_df(sim_params, datapath, savepath, worker_id=worker_id)
    return gnss_const, signal_family, df_outs


if __name__ == "__main__":
    results = {}

    # adjust max_workers as you like (None -> use number of CPUs)
    with ProcessPoolExecutor(max_workers=None) as executor:
        # submit all jobs
        future_to_combo = {
            executor.submit(run_single, combo, worker_id=i): combo for i, combo in enumerate(COMBOS)
        }

        for future in as_completed(future_to_combo):
            gnss_const, signal_family, df_outs = future.result()
            results[(gnss_const, signal_family)] = df_outs
            print(f"Finished {gnss_const}, signal family {signal_family}")
