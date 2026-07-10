#!/usr/bin/env python
"""Run the agent-based lunar ephemeris/almanac datasize-accuracy study and export to .npz.

The datasize/accuracy sweep (propagate a truth arc, fit + quantize the broadcast models
over several fitting-window lengths, size the bit budget) is compute-heavy; running it here
rather than in the notebook keeps the tutorial fast. This script builds a ``pnt.Simulation``
from ``configs/ephemeris.yaml`` -- a thin ``EphemerisManager`` agent hosting an
``EphemerisApp`` -- and runs it twice on the same ELFO orbit:

    eph : short (<= 2 h) windows on a 1-day, 60-s truth arc   (keep CartesianEphemeris)
    alm : long (<= 15 d) windows on a 30-day, 600-s truth arc (keep Almanac)

Each run's results are read from the coordinator app
(``sim.get_agent("EphemerisManager").get_application().get_results()``) and written to one
``.npz`` per sweep (truth arc in MOON_PA + per-window results), plus ``meta.json`` that
``ex9_ephemeris.ipynb`` reads back for plotting.

Usage:
    python ex9_run_ephemeris.py [--config PATH] [--outdir DIR]
"""
import argparse
import copy
import json
import os
import sys
from pathlib import Path

import numpy as np
import yaml

SECS_DAY = 86400.0
SECS_MINUTE = 60.0
MIN_PER_DAY = SECS_DAY / SECS_MINUTE  # 1440.0


def find_repo_root() -> Path:
    for b in [Path.cwd(), *Path.cwd().parents]:
        if (b / "python/pylupnt/__init__.py").exists():
            return b
    raise RuntimeError("Could not locate repo root (python/pylupnt/__init__.py)")


def window_results_to_dict(results) -> dict:
    """Stack a list[EphemerisWindowResult] into column arrays."""
    return {
        "fit_window_min": np.array([r.fit_window_min for r in results]),
        "num_params": np.array([r.num_params for r in results]),
        "total_bits": np.array([r.total_bits for r in results]),
        "pos_rms_m": np.array([r.pos_rms_m for r in results]),
        "pos_p95_m": np.array([r.pos_p95_m for r in results]),
        "vel_rms_mps": np.array([r.vel_rms_mps for r in results]),
        "vel_p95_mps": np.array([r.vel_p95_mps for r in results]),
    }


def main() -> None:
    repo = find_repo_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(repo / "configs" / "ephemeris.yaml"))
    ap.add_argument("--outdir", default=str(repo / "output" / "python_examples" / "ex9_data"))
    args = ap.parse_args()

    sys.path.insert(0, str((repo / "python").resolve()))
    if (repo / "data/LuPNT_data").is_dir():
        os.environ.setdefault("LUPNT_DATA_PATH", str((repo / "data/LuPNT_data").resolve()))
    import pylupnt as pnt  # noqa: E402

    with open(args.config) as f:
        base_cfg = yaml.safe_load(f)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Two sweeps on the same ELFO orbit: the ephemeris over short windows (base config)
    # and the almanac over long windows (30-day, coarser arc). Each keeps only its own
    # model's results (see the notebook).
    variants = {
        "eph": {},  # base config = ephemeris sweep
        "alm": {
            "duration_days": 30.0,
            "sample_dt_s": 600.0,
            "fit_window_minutes": [1.0 * MIN_PER_DAY, 7.0 * MIN_PER_DAY, 15.0 * MIN_PER_DAY],
        },
    }

    meta = {"C_mps": float(pnt.C), "R_MOON_m": float(pnt.R_MOON), "GM_MOON": float(pnt.GM_MOON)}
    for tag, overrides in variants.items():
        cfg = copy.deepcopy(base_cfg)
        acfg = cfg["agents"]["EphemerisManager"]["application"]
        acfg.update(overrides)
        print(f"Running ephemeris study ({tag}): "
              f"duration_days={acfg['duration_days']}, windows={acfg['fit_window_minutes']}")
        sim = pnt.Simulation(cfg)
        sim.run()
        app = sim.get_agent("EphemerisManager").get_application()
        res = app.get_results()

        d = {
            "t_truth_s": np.asarray(res.t_truth_s).reshape(-1),
            "rv_truth": np.asarray(res.rv_truth),  # [N x 6], MOON_PA (output_frame)
        }
        for k, v in window_results_to_dict(res.cartesian_results).items():
            d[f"cart_{k}"] = v
        for k, v in window_results_to_dict(res.almanac_results).items():
            d[f"alm_{k}"] = v
        np.savez_compressed(outdir / f"results_{tag}.npz", **d)

        cart = res.cartesian_results
        alm = res.almanac_results
        print(f"  {tag}: cart p95 pos [m] = {[round(r.pos_p95_m, 4) for r in cart]}, "
              f"alm p95 pos [m] = {[round(r.pos_p95_m, 1) for r in alm]}")

        if tag == "eph":
            oc = acfg["orbit"]
            meta.update({
                "start_epoch_utc": acfg["start_epoch_utc"],
                "orbit": oc,
                "output_frame": acfg["output_frame"],
                "cartesian_poly_order": acfg["cartesian_poly_order"],
                "almanac_poly_order": acfg["almanac_poly_order"],
                "almanac_num_fourier_terms": acfg["almanac_num_fourier_terms"],
                "datasize_precision_m": acfg["datasize_precision_m"],
                "num_windows": acfg["num_windows"],
                "eph_sample_dt_s": acfg["sample_dt_s"],
            })
        if tag == "alm":
            meta["alm_sample_dt_s"] = acfg["sample_dt_s"]

    with open(outdir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote results to {outdir}")


if __name__ == "__main__":
    main()
