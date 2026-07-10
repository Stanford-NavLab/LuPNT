#!/usr/bin/env python
"""Run the agent-based distributed ISL ODTS simulation and export results to .npz.

The five parallel Schmidt-EKF filters (plus the centralized ground filter) over the
6-hour arc are compute-heavy; running them here (rather than inside the notebook) keeps
the tutorial notebook fast. This script builds a ``pnt.Simulation`` from
``configs/isl_odts.yaml`` and runs three variants on the same truth/geometry:

    full     : station beacon network + 10-min inter-agent (consider-state) exchange
    noex     : exchange OFF  (each filter keeps its static consider states)
    nost     : stations OFF  (crosslinks + exchange only, no absolute timing anchor)

It reads each run's results from the coordinator app
(``sim.get_agent("IslManager").get_application().get_results()``) and writes one ``.npz``
per variant (all IslOdtsResults series) plus ``meta.json`` that
``ex8_isl_odts.ipynb`` reads back for plotting.

Usage:
    python ex8_run_isl_odts.py [--config PATH] [--outdir DIR] [--duration HH:MM:SS]
"""
import argparse
import copy
import json
import os
import sys
from pathlib import Path

import numpy as np
import yaml


def find_repo_root() -> Path:
    for b in [Path.cwd(), *Path.cwd().parents]:
        if (b / "python/pylupnt/__init__.py").exists():
            return b
    raise RuntimeError("Could not locate repo root (python/pylupnt/__init__.py)")


def _hms_to_s(s: str) -> float:
    parts = [float(p) for p in s.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


def results_to_dict(res) -> dict:
    """Flatten an IslOdtsResults into a dict of numpy arrays (stacking the per-satellite
    / per-station lists along a leading axis) for np.savez."""
    n_sat = len(res.satellite_names)

    def stack(lst):
        return np.stack([np.asarray(a) for a in lst], axis=0) if len(lst) else np.zeros((0,))

    d = {
        "t_s": np.asarray(res.t_s).reshape(-1),
        "truth_states": stack(res.truth_states),          # [n_sat, N, 8]
        "est": stack(res.est),                            # [n_sat, N, 8*n_sat]
        "cov_diag": stack(res.cov_diag),                  # [n_sat, N, 8*n_sat]
        "cov_own_full": stack(res.cov_own_full),          # [n_sat, N, 64]
        "range_true_m": np.asarray(res.range_true_m),
        "range_rate_true_mps": np.asarray(res.range_rate_true_mps),
        "range_obs_m": np.asarray(res.range_obs_m),
        "range_rate_obs_mps": np.asarray(res.range_rate_obs_mps),
        "range_resid_m": stack(res.range_resid_m),        # [n_sat, N, n_links]
        "cn0_dbhz": np.asarray(res.cn0_dbhz),
        "time_transfer_true_m": np.asarray(res.time_transfer_true_m),
        "time_transfer_obs_m": np.asarray(res.time_transfer_obs_m),
        "station_pos_mci": stack(res.station_pos_mci),    # [n_stn, N, 3]
        "station_visible": np.asarray(res.station_visible),
        "station_pr_true_m": np.asarray(res.station_pr_true_m),
        "station_pr_obs_m": np.asarray(res.station_pr_obs_m),
        "station_pr_resid_m": np.asarray(res.station_pr_resid_m),
        "est_central": np.asarray(res.est_central),
        "cov_central_full": stack(res.cov_central_full),  # [n_sat, N, 64]
        "n_sat": np.asarray(n_sat),
    }
    return d


def main() -> None:
    repo = find_repo_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(repo / "configs" / "isl_odts.yaml"))
    ap.add_argument("--outdir", default=str(repo / "output" / "python_examples" / "ex8_data"))
    ap.add_argument("--duration", default=None, help="Override scenario duration, e.g. 06:00:00")
    args = ap.parse_args()

    sys.path.insert(0, str((repo / "python").resolve()))
    if (repo / "data/LuPNT_data").is_dir():
        os.environ.setdefault("LUPNT_DATA_PATH", str((repo / "data/LuPNT_data").resolve()))
    import pylupnt as pnt  # noqa: E402

    with open(args.config) as f:
        base_cfg = yaml.safe_load(f)
    app_path = base_cfg["agents"]["IslManager"]["application"]
    if args.duration:
        base_cfg["duration"] = args.duration
        app_path["duration_s"] = _hms_to_s(args.duration)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = {
        "full": {},
        "noex": {"consider_exchange_interval_s": 0.0},
        "nost": {"surface_stations": [], "enable_centralized_ground_filter": False},
        # Consider-state exchange via Covariance Intersection (consistent fusion) with a
        # matched (larger) process noise, instead of the naive-overwrite baseline -- the
        # Section 7-8 consistency (NEES) comparison.
        "ci": {"exchange_use_covariance_intersection": True, "exchange_ci_weight": 0.95,
               "process_accel_sigma_mps2": 1.0e-3},
        # Crosslink clock-observable sweep (CI exchange, weight 0.9): range+range-rate only,
        # then + two-way time transfer, then + frequency transfer -- Section 9.
        "tf_rr": {"exchange_use_covariance_intersection": True, "exchange_ci_weight": 0.9,
                  "enable_two_way_time_transfer": False, "enable_two_way_frequency_transfer": False},
        "tf_tt": {"exchange_use_covariance_intersection": True, "exchange_ci_weight": 0.9,
                  "enable_two_way_time_transfer": True, "enable_two_way_frequency_transfer": False},
        "tf_ttft": {"exchange_use_covariance_intersection": True, "exchange_ci_weight": 0.9,
                    "enable_two_way_time_transfer": True, "enable_two_way_frequency_transfer": True},
    }

    meta = {}
    for tag, overrides in variants.items():
        cfg = copy.deepcopy(base_cfg)
        acfg = cfg["agents"]["IslManager"]["application"]
        acfg.update(overrides)
        print(f"Running distributed ISL ODTS ({tag}): duration={cfg['duration']}")
        sim = pnt.Simulation(cfg)
        sim.run()
        app = sim.get_agent("IslManager").get_application()
        res = app.get_results()

        d = results_to_dict(res)
        np.savez_compressed(outdir / f"results_{tag}.npz", **d)

        if tag == "full":
            n_sat = len(res.satellite_names)
            meta = {
                "duration": cfg["duration"],
                "dt_s": float(acfg["dt_s"]),
                "start_epoch_utc": acfg["start_epoch_utc"],
                "satellite_names": list(res.satellite_names),
                "n_sat": n_sat,
                "n_links": n_sat - 1,
                "consider_exchange_interval_s": float(acfg["consider_exchange_interval_s"]),
                "station_names": [s.get("name", f"GS{i}")
                                  for i, s in enumerate(acfg.get("surface_stations", []))],
                "C_mps": float(pnt.C),
            }
        # Console summary: final own position error per satellite.
        final_pos = []
        for j in range(len(res.satellite_names)):
            own = np.asarray(res.est[j])[-1, 8 * j:8 * j + 3]
            tru = np.asarray(res.truth_states[j])[-1, :3]
            final_pos.append(float(np.linalg.norm(own - tru)))
        meta[f"final_pos_err_{tag}"] = final_pos
        print(f"  mean final own pos err ({tag}): {np.mean(final_pos):.1f} m")

    with open(outdir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote results to {outdir}")


if __name__ == "__main__":
    main()
