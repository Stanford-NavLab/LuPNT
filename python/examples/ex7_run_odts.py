#!/usr/bin/env python
"""Run the agent-based ground-station ODTS simulation and export results to CSV.

The orbit-determination arc (36 h by default) runs the centralized batch + SRIF
filter, which is compute-heavy. Running it here (rather than inside the notebook)
keeps the tutorial notebook fast and robust: this script does the work once and
writes plain CSV files that ``ex7_groundstation_odts.ipynb`` reads back for
plotting.

Usage:
    python ex7_run_odts.py [--config PATH] [--outdir DIR] [--duration HH:MM:SS]

Outputs (under --outdir, default output/python_examples/ex7_data/):
    grid.csv          uniform epoch grid: truth/estimate/SRIF states + covariances
    iterations.csv    batch-filter iteration history
    elevation.csv     per-station topocentric elevation series (long format)
    measurements.csv  per-station range/range-rate observations (long format)
    epoch_solution.csv  true / initial-guess / estimated epoch state
    covariance.csv    6x6 formal covariance of the estimated epoch state
    meta.json         scalars + config (converged, iterations, noise, labels, ...)
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def find_repo_root() -> Path:
    for b in [Path.cwd(), *Path.cwd().parents]:
        if (b / "python/pylupnt/__init__.py").exists():
            return b
    raise RuntimeError("Could not locate repo root (python/pylupnt/__init__.py)")


def main() -> None:
    repo = find_repo_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(repo / "configs" / "ground_station_odts.yaml"))
    ap.add_argument("--outdir", default=str(repo / "output" / "python_examples" / "ex7_data"))
    ap.add_argument(
        "--duration",
        default=None,
        help="Override the scenario duration, e.g. 12:00:00 (default: config value)",
    )
    args = ap.parse_args()

    sys.path.insert(0, str((repo / "python").resolve()))
    if (repo / "data/LuPNT_data/ephemeris").is_dir():
        os.environ.setdefault("LUPNT_DATA_PATH", str((repo / "data/LuPNT_data").resolve()))
    import pylupnt as pnt  # noqa: E402  (after sys.path/env setup)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    if args.duration:
        cfg["duration"] = args.duration

    station_order = [n for n in cfg["agents"] if cfg["agents"][n]["class"] == "GroundStation"]
    track_cfg = cfg["agents"][station_order[0]]["application"]

    print(f"Running ground-station ODTS: duration={cfg['duration']}, " f"stations={station_order}")
    sim = pnt.Simulation(cfg)
    sim.run()

    mgr = sim.get_agent("gs_manager").get_application()
    stations = {n: sim.get_agent(n).get_application() for n in station_order}

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Grid time series ----
    t_grid = np.asarray(mgr.time_grid()).reshape(-1)
    N = len(t_grid)
    comp = ["rx", "ry", "rz", "vx", "vy", "vz"]
    cov_cols = [f"{r}{c}" for r in range(6) for c in range(6)]

    grid = {"t_s": t_grid}
    for pre, arr in [
        ("tr", mgr.truth_state()),
        ("es", mgr.estimated_state()),
        ("sf", mgr.srif_filtered_state()),
        ("ss", mgr.srif_smoothed_state()),
    ]:
        A = np.asarray(arr).reshape(N, 6)
        for j, c in enumerate(comp):
            grid[f"{pre}_{c}"] = A[:, j]
    for pre, arr in [
        ("ec", mgr.estimated_covariance()),
        ("sfc", mgr.srif_filtered_covariance()),
        ("ssc", mgr.srif_smoothed_covariance()),
    ]:
        A = np.asarray(arr).reshape(N, 36)
        for j, c in enumerate(cov_cols):
            grid[f"{pre}_{c}"] = A[:, j]
    pd.DataFrame(grid).to_csv(outdir / "grid.csv", index=False)

    # ---- Batch iteration history ----
    pd.DataFrame(
        {
            "iter": np.arange(mgr.num_iterations()),
            "pos_error_m": np.asarray(mgr.iteration_pos_error()).reshape(-1),
            "vel_error_mps": np.asarray(mgr.iteration_vel_error()).reshape(-1),
            "correction_norm": np.asarray(mgr.iteration_correction_norm()).reshape(-1),
            "weighted_rms": np.asarray(mgr.iteration_weighted_rms()).reshape(-1),
            "rms_range_m": np.asarray(mgr.iteration_rms_range()).reshape(-1),
            "rms_range_rate_mps": np.asarray(mgr.iteration_rms_range_rate()).reshape(-1),
        }
    ).to_csv(outdir / "iterations.csv", index=False)

    # ---- Per-station observation series (long format) ----
    elev_rows, meas_rows = [], []
    for n in station_order:
        app = stations[n]
        et, ed = np.asarray(app.elevation_time()), np.asarray(app.elevation_deg())
        for t, e in zip(et, ed):
            elev_rows.append({"station": n, "t_s": float(t), "elevation_deg": float(e)})
        mt = np.asarray(app.measurement_time())
        mr, mrr = np.asarray(app.measurement_range()), np.asarray(app.measurement_range_rate())
        for t, r, rr in zip(mt, mr, mrr):
            meas_rows.append(
                {"station": n, "t_s": float(t), "range_m": float(r), "range_rate_mps": float(rr)}
            )
    pd.DataFrame(elev_rows).to_csv(outdir / "elevation.csv", index=False)
    pd.DataFrame(meas_rows).to_csv(outdir / "measurements.csv", index=False)

    # ---- Epoch solution + covariance ----
    pd.DataFrame(
        {
            "which": ["true", "initial_guess", "estimated"],
            **{
                c: [
                    np.asarray(mgr.x0_true()).reshape(-1)[j],
                    np.asarray(mgr.x0_initial_guess()).reshape(-1)[j],
                    np.asarray(mgr.x0_estimated()).reshape(-1)[j],
                ]
                for j, c in enumerate(comp)
            },
        }
    ).to_csv(outdir / "epoch_solution.csv", index=False)
    pd.DataFrame(np.asarray(mgr.covariance()).reshape(6, 6), columns=comp, index=comp).to_csv(
        outdir / "covariance.csv"
    )

    # ---- Scalars + config metadata ----
    labels = {"DSS14": "DSS14 (Goldstone)", "DSS43": "DSS43 (Canberra)", "DSS63": "DSS63 (Madrid)"}
    obs_interval_s = float(cfg["agents"]["gs_manager"]["application"].get("obs_interval_s", 300.0))
    meta = {
        "duration": cfg["duration"],
        "obs_interval_s": obs_interval_s,
        "converged": bool(mgr.converged()),
        "num_iterations": int(mgr.num_iterations()),
        "num_measurements": int(mgr.num_measurements()),
        "station_order": station_order,
        "station_labels": {n: labels.get(n, n) for n in station_order},
        "elevation_mask_deg": float(track_cfg.get("elevation_mask_deg", 10.0)),
        "range_sigma_m": float(track_cfg.get("range_sigma_m", 10.0)),
        "range_rate_sigma_mps": float(track_cfg.get("range_rate_sigma_mps", 1.0e-3)),
    }
    with open(outdir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote CSVs to {outdir}")
    print(
        f"  converged={meta['converged']} iters={meta['num_iterations']} "
        f"measurements={meta['num_measurements']}"
    )
    pe = np.linalg.norm(
        np.asarray(mgr.x0_estimated()).reshape(-1)[:3] - np.asarray(mgr.x0_true()).reshape(-1)[:3]
    )
    print(f"  final position error: {pe:.3f} m")


if __name__ == "__main__":
    main()
