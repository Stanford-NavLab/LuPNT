#!/usr/bin/env python
"""Run the physical-agent distributed ISL ODTS scenario and export results to .npz.

Each satellite is a physical ``Spacecraft`` agent hosting a ``SatelliteOdtsApp`` onboard
Schmidt-EKF; each surface station is a physical ``SurfaceStation`` running a
``StationBeaconSensor`` that feeds a ``SurfaceStationManager``'s ``GroundOdtsApp`` (the
centralized ground filter, no ISL). This script builds a ``pnt.Simulation`` from
``configs/isl_odts_distributed.yaml`` and runs several variants on the same truth/geometry:

    full     : station beacon network + inter-agent (consider-state) exchange
    noex     : exchange OFF  (each filter keeps its static consider states)
    nost     : stations OFF  (crosslinks + exchange only, no absolute timing anchor)
    ci       : Covariance-Intersection exchange (consistent fusion) + matched process noise
    tf_rr / tf_tt / tf_ttft : crosslink clock-observable sweep (range/rate; +time; +freq transfer)

Per-satellite onboard results are read from each ``SatelliteOdtsApp`` (``own_estimate`` /
``truth_state`` / ``own_cov_full``) and the ground filter from the ``GroundOdtsApp``
(``est_central`` / ``cov_central_full``); one ``.npz`` per variant + ``meta.json`` are written
for ``ex8_isl_odts.ipynb``.

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


def sat_names(cfg) -> list:
    return [
        n
        for n, a in cfg["agents"].items()
        if a.get("application", {}).get("class") == "SatelliteOdtsApp"
    ]


def apply_variant(cfg, sat_over: dict, drop_stations: bool = False):
    """Apply per-satellite app overrides to every SatelliteOdtsApp; optionally remove the
    surface-station agents + their ground manager (the `nost` variant)."""
    for name in list(cfg["agents"]):
        agent = cfg["agents"][name]
        app = agent.get("application", {})
        if app.get("class") == "SatelliteOdtsApp":
            app.update(sat_over)
            if drop_stations:
                app["stations"] = []
        if drop_stations and app.get("class") in ("StationBeaconSensor", "GroundOdtsApp"):
            del cfg["agents"][name]


def extract(sim, names) -> dict:
    apps = [sim.get_agent(sn).get_application() for sn in names]
    d = {
        "t_s": np.asarray(apps[0].time_grid()).reshape(-1),
        "truth_states": np.stack([np.asarray(a.truth_state()) for a in apps]),  # [n_sat, N, 8]
        "own_est": np.stack([np.asarray(a.own_estimate()) for a in apps]),  # [n_sat, N, 8]
        "own_cov_full": np.stack([np.asarray(a.own_cov_full()) for a in apps]),  # [n_sat, N, 64]
        "n_sat": np.asarray(len(names)),
    }
    try:
        g = sim.get_agent("gs_manager").get_application()
        d["est_central"] = np.asarray(g.est_central())  # [N, 8*n_sat]
        d["cov_central_full"] = np.stack([np.asarray(c) for c in g.cov_central_full()])
    except Exception:
        pass  # `nost` variant has no ground manager
    return d


def main() -> None:
    repo = find_repo_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(repo / "configs" / "isl_odts_distributed.yaml"))
    ap.add_argument("--outdir", default=str(repo / "output" / "python_examples" / "ex8_data"))
    ap.add_argument("--duration", default=None, help="Override scenario duration, e.g. 06:00:00")
    args = ap.parse_args()

    sys.path.insert(0, str((repo / "python").resolve()))
    if (repo / "data/LuPNT_data").is_dir():
        os.environ.setdefault("LUPNT_DATA_PATH", str((repo / "data/LuPNT_data").resolve()))
    import pylupnt as pnt  # noqa: E402

    with open(args.config) as f:
        base_cfg = yaml.safe_load(f)
    dur_s = None
    if args.duration:
        base_cfg["duration"] = args.duration
        dur_s = _hms_to_s(args.duration)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ci = {"exchange_use_covariance_intersection": True, "exchange_ci_weight": 0.9}
    variants = {
        "full": ({}, False),
        "noex": ({"consider_exchange_interval_s": 0.0}, False),
        "nost": ({}, True),
        "ci": (
            {
                "exchange_use_covariance_intersection": True,
                "exchange_ci_weight": 0.95,
                "process_accel_sigma_mps2": 1.0e-3,
            },
            False,
        ),
        "tf_rr": (
            {
                **ci,
                "enable_two_way_time_transfer": False,
                "enable_two_way_frequency_transfer": False,
            },
            False,
        ),
        "tf_tt": (
            {
                **ci,
                "enable_two_way_time_transfer": True,
                "enable_two_way_frequency_transfer": False,
            },
            False,
        ),
        "tf_ttft": (
            {**ci, "enable_two_way_time_transfer": True, "enable_two_way_frequency_transfer": True},
            False,
        ),
    }

    names = sat_names(base_cfg)
    meta = {
        "satellite_names": names,
        "n_sat": len(names),
        "n_links": len(names) - 1,
        "duration": base_cfg["duration"],
        "C_mps": float(pnt.C),
    }
    for tag, (sat_over, drop_st) in variants.items():
        cfg = copy.deepcopy(base_cfg)
        if dur_s is not None:
            for a in cfg["agents"].values():
                if a.get("application", {}).get("duration_s"):
                    a["application"]["duration_s"] = dur_s
        apply_variant(cfg, sat_over, drop_st)
        print(f"Running distributed ISL ODTS ({tag}): duration={cfg['duration']}")
        sim = pnt.Simulation(cfg)
        sim.run()
        d = extract(sim, names)
        np.savez_compressed(outdir / f"results_{tag}.npz", **d)

        # Console summary: final onboard own position error per satellite.
        final = [
            float(np.linalg.norm(d["own_est"][j, -1, :3] - d["truth_states"][j, -1, :3]))
            for j in range(len(names))
        ]
        meta[f"final_pos_err_{tag}"] = final
        print(f"  mean final own pos err ({tag}): {np.mean(final):.1f} m")

    with open(outdir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote results to {outdir}")


if __name__ == "__main__":
    main()
