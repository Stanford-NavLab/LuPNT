#!/usr/bin/env python
"""Run the agent-based surface-rover navigation simulation and export results to CSV.

The strapdown-INS EKF over the ~30-min driving arc is compute-heavy; running it here
(rather than inside the notebook) keeps the tutorial notebook fast and robust. This
script builds a ``pnt.Simulation`` from ``configs/surface_rover_nav.yaml``, runs it
with the DEM altitude constraint enabled and (as an ablation) disabled, and writes
plain CSV / NumPy files that ``ex10_surface_rover.ipynb`` reads back for plotting.

Usage:
    python ex10_run_surface_nav.py [--config PATH] [--outdir DIR] [--duration HH:MM:SS]

Outputs (under --outdir, default output/python_examples/ex10_data/):
    series_dem_on.csv    per-epoch truth/estimate error + covariance (DEM constraint on)
    series_dem_off.csv   same, with the DEM altitude constraint disabled (ablation)
    dem_x.npy dem_y.npy dem_elevation.npy   terrain grid (native projected meters)
    meta.json            site info, satellite names, scalars, config
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


def series_frame(app) -> pd.DataFrame:
    """Flatten a SurfaceRoverNavApp's result accessors into a tidy per-epoch DataFrame."""
    d = {"t_s": np.asarray(app.time_s()).reshape(-1)}
    v3 = {
        "pos_err": app.pos_err_enu(), "pos_sigma": app.pos_sigma_enu(),
        "accel_bias_err": app.accel_bias_err(), "accel_bias_sigma": app.accel_bias_sigma(),
        "gyro_bias_err": app.gyro_bias_err(), "gyro_bias_sigma": app.gyro_bias_sigma(),
        "att_err_deg": app.att_err_deg(), "att_sigma_deg": app.att_sigma_deg(),
    }
    for name, arr in v3.items():
        A = np.asarray(arr).reshape(-1, 3)
        for j, c in enumerate(["x", "y", "z"]):
            d[f"{name}_{c}"] = A[:, j]
    for name, arr in [("pos_err_norm", app.pos_err_norm()),
                      ("clock_bias_err", app.clock_bias_err()),
                      ("clock_bias_sigma", app.clock_bias_sigma()),
                      ("rover_alt_truth", app.rover_alt_truth())]:
        d[name] = np.asarray(arr).reshape(-1)
    d["n_visible"] = np.asarray(app.n_visible()).reshape(-1)
    for name, arr in [("track_truth", app.rover_track_enu_truth()),
                      ("track_est", app.rover_track_enu_est())]:
        A = np.asarray(arr).reshape(-1, 2)
        d[f"{name}_e"], d[f"{name}_n"] = A[:, 0], A[:, 1]
    return pd.DataFrame(d)


def main() -> None:
    repo = find_repo_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(repo / "configs" / "surface_rover_nav.yaml"))
    ap.add_argument("--outdir", default=str(repo / "output" / "python_examples" / "ex10_data"))
    ap.add_argument("--duration", default=None, help="Override scenario duration, e.g. 00:15:00")
    args = ap.parse_args()

    sys.path.insert(0, str((repo / "python").resolve()))
    if (repo / "data/LuPNT_data").is_dir():
        os.environ.setdefault("LUPNT_DATA_PATH", str((repo / "data/LuPNT_data").resolve()))
    os.environ.setdefault("LUPNT_SKIP_DEM_DOWNLOAD", "1")
    import pylupnt as pnt  # noqa: E402

    with open(args.config) as f:
        base_cfg = yaml.safe_load(f)
    if args.duration:
        base_cfg["duration"] = args.duration
        base_cfg["agents"]["Rover"]["application"]["duration_s"] = _hms_to_s(args.duration)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dem_saved = False
    meta = {}
    for tag, dem_on in [("dem_on", True), ("dem_off", False)]:
        cfg = yaml.safe_load(yaml.safe_dump(base_cfg))  # deep copy
        cfg["agents"]["Rover"]["application"]["enable_dem_constraint"] = dem_on
        print(f"Running surface-rover nav ({tag}): duration={cfg['duration']}")
        sim = pnt.Simulation(cfg)
        sim.run()
        app = sim.get_agent("Rover").get_application()

        series_frame(app).to_csv(outdir / f"series_{tag}.csv", index=False)
        if not dem_saved:
            np.save(outdir / "dem_x.npy", np.asarray(app.dem_x()))
            np.save(outdir / "dem_y.npy", np.asarray(app.dem_y()))
            np.save(outdir / "dem_elevation.npy", np.asarray(app.dem_elevation()))
            acfg = cfg["agents"]["Rover"]["application"]
            meta = {
                "duration": cfg["duration"],
                "dt_s": float(acfg["dt_s"]),
                "site_id": app.site_id(),
                "site_name": app.site_name(),
                "site_lat_deg": float(cfg["world"]["dem"]["site_lat_deg"]),
                "site_lon_deg": float(cfg["world"]["dem"]["site_lon_deg"]),
                "satellite_names": list(app.satellite_names()),
                "pseudorange_sigma_m": float(acfg["pseudorange_sigma_m"]),
                "sise_m": float(acfg["sise_m"]),
                "elevation_mask_deg": float(acfg["elevation_mask_deg"]),
                "dem_sigma_m": float(acfg["dem_sigma_m"]),
            }
            dem_saved = True
        pe = float(np.asarray(app.pos_err_norm()).reshape(-1)[-1])
        up = float(abs(np.asarray(app.pos_err_enu()).reshape(-1, 3)[-1, 2]))
        meta[f"final_pos_err_norm_{tag}"] = pe
        meta[f"final_up_err_{tag}"] = up
        print(f"  final 3D position error ({tag}): {pe:.3f} m  (Up error {up:.3f} m)")

    with open(outdir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote CSVs to {outdir}")


def _hms_to_s(s: str) -> float:
    parts = [float(p) for p in s.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


if __name__ == "__main__":
    main()
