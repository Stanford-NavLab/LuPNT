#!/usr/bin/env python
"""Run the agent-based lunar-lander descent-navigation simulation and export results to CSV.

Builds a ``pnt.Simulation`` from ``configs/lander_nav.yaml`` (a thin ``Lander`` agent
hosting a ``LanderNavApp`` MEKF that fuses IMU + radar altimeter + crater bearings +
LunaNet pseudoranges), runs the full-sensor descent plus a sensor ablation (each aiding
source disabled in turn), and writes CSV / NumPy files that ``ex11_lander_navigation.ipynb``
reads back for plotting.

Usage:
    python ex11_run_lander_nav.py [--config PATH] [--outdir DIR] [--duration HH:MM:SS]

Outputs (under --outdir, default output/python_examples/ex11_data/):
    series_all.csv       per-epoch truth/estimate error + covariance (all sensors)
    craters.csv          synthetic crater map (East, North) [m]
    dem_x.npy dem_y.npy dem_elevation.npy   terrain grid (native projected meters)
    ablation.csv         final 3D position error with each aiding sensor disabled
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
    d = {"t_s": np.asarray(app.time_s()).reshape(-1)}
    v3 = {
        "pos_err": app.pos_err_enu(), "pos_sigma": app.pos_sigma_enu(),
        "accel_bias_err": app.accel_bias_err(), "accel_bias_sigma": app.accel_bias_sigma(),
        "gyro_bias_err": app.gyro_bias_err(), "gyro_bias_sigma": app.gyro_bias_sigma(),
        "att_err_deg": app.att_err_deg(), "att_sigma_deg": app.att_sigma_deg(),
        "traj_truth": app.traj_enu_truth(), "traj_est": app.traj_enu_est(),
    }
    for name, arr in v3.items():
        A = np.asarray(arr).reshape(-1, 3)
        for j, c in enumerate(["x", "y", "z"]):
            d[f"{name}_{c}"] = A[:, j]
    for name, arr in [("pos_err_norm", app.pos_err_norm()), ("vel_err_norm", app.vel_err_norm()),
                      ("clock_bias_err", app.clock_bias_err()),
                      ("clock_bias_sigma", app.clock_bias_sigma()),
                      ("alt_truth", app.alt_truth()), ("alt_est", app.alt_est())]:
        d[name] = np.asarray(arr).reshape(-1)
    d["n_visible_sat"] = np.asarray(app.n_visible_sat()).reshape(-1)
    d["n_craters"] = np.asarray(app.n_craters()).reshape(-1)
    return pd.DataFrame(d)


def run(cfg, pnt):
    sim = pnt.Simulation(yaml.safe_load(yaml.safe_dump(cfg)))
    sim.run()
    return sim.get_agent("Lander").get_application()


def main() -> None:
    repo = find_repo_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(repo / "configs" / "lander_nav.yaml"))
    ap.add_argument("--outdir", default=str(repo / "output" / "python_examples" / "ex11_data"))
    ap.add_argument("--duration", default=None, help="Override scenario duration, e.g. 00:05:00")
    args = ap.parse_args()

    sys.path.insert(0, str((repo / "python").resolve()))
    if (repo / "data/LuPNT_data").is_dir():
        os.environ.setdefault("LUPNT_DATA_PATH", str((repo / "data/LuPNT_data").resolve()))
    os.environ.setdefault("LUPNT_SKIP_DEM_DOWNLOAD", "1")
    import pylupnt as pnt  # noqa: E402

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    if args.duration:
        cfg["duration"] = args.duration
        cfg["agents"]["Lander"]["application"]["duration_s"] = _hms_to_s(args.duration)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    acfg = cfg["agents"]["Lander"]["application"]

    print(f"Running lander nav (all sensors): duration={cfg['duration']}")
    app = run(cfg, pnt)
    series_frame(app).to_csv(outdir / "series_all.csv", index=False)
    craters = np.asarray(app.crater_enu()).reshape(-1, 2)
    pd.DataFrame({"east_m": craters[:, 0], "north_m": craters[:, 1]}).to_csv(
        outdir / "craters.csv", index=False)
    np.save(outdir / "dem_x.npy", np.asarray(app.dem_x()))
    np.save(outdir / "dem_y.npy", np.asarray(app.dem_y()))
    np.save(outdir / "dem_elevation.npy", np.asarray(app.dem_elevation()))
    pe_all = float(np.asarray(app.pos_err_norm()).reshape(-1)[-1])

    # ---- Sensor ablation: disable one aiding source at a time ----
    abl = [("all sensors", {})]
    abl += [("no craters", {"enable_craters": False}),
            ("no altimeter", {"enable_altimeter": False}),
            ("no LunaNet", {"enable_lunanet": False})]
    rows = []
    for label, over in abl:
        c = yaml.safe_load(yaml.safe_dump(cfg))
        c["agents"]["Lander"]["application"].update(over)
        a = run(c, pnt)
        pe = float(np.asarray(a.pos_err_norm()).reshape(-1)[-1])
        rows.append({"config": label, "final_pos_err_norm_m": pe})
        print(f"  {label:14s}: final 3D position error {pe:.3f} m")
    pd.DataFrame(rows).to_csv(outdir / "ablation.csv", index=False)

    meta = {
        "duration": cfg["duration"],
        "dt_s": float(acfg["dt_s"]),
        "site_id": app.site_id(),
        "site_name": app.site_name(),
        "site_lat_deg": float(cfg["world"]["dem"]["site_lat_deg"]),
        "site_lon_deg": float(cfg["world"]["dem"]["site_lon_deg"]),
        "satellite_names": list(app.satellite_names()),
        "altimeter_sigma_m": float(acfg["altimeter_sigma_m"]),
        "crater_sigma_arcsec": float(acfg["crater_sigma_arcsec"]),
        "pseudorange_sigma_m": float(acfg["pseudorange_sigma_m"]),
        "final_pos_err_norm_all": pe_all,
        "touchdown_alt_truth_m": float(np.asarray(app.alt_truth()).reshape(-1)[-1]),
    }
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
