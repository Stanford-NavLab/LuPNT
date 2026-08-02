"""Run the ex6 lunar GNSS ODTS EKF via the agent-based ``pnt.Simulation`` and export results.

This is the standalone equivalent of the notebook cell that runs the EKF in
``ex6_gnss_odts.ipynb``. It drives the agent-based architecture: a physical ``receiver`` Spacecraft
hosts a ``LunarGnssOdtsApp`` whose scheduled step runs the whole Monte-Carlo
ODTS body (numerically bit-identical to the former ``LunarGnssODTSSimulation``).

Run the slower precompute stages first (they build the Stage 1 link cache + Stage 2 plasma
delays consumed here)::

    pixi run python python/examples/ex6_precompute.py
    pixi run python python/examples/ex6_run_gnss_odts.py

The EKF writes ``trajectory_mc<N>.csv`` + ``summary.csv`` under the output dir (the C++ engine
does this); this script additionally writes ``ex6_results.npz`` (per-seed summaries) for quick
programmatic access. The notebook reads the trajectory CSV to plot.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def log(message=""):
    print(message, flush=True)


def _require_input_file(path, label):
    path = Path(path)
    if path.exists() and path.stat().st_size > 0:
        return
    raise FileNotFoundError(
        f"Missing {label}: {path}\n"
        "Run the precompute stages first:\n"
        "  pixi run python python/examples/ex6_precompute.py"
    )


def _print_summary(summary, app_cfg):
    log(f"Final position error:    {summary.final_position_error_m:8.3f} m")
    log(f"Final velocity error:    {summary.final_velocity_error_mps * 1e3:8.3f} mm/s")
    log(f"Final clock bias error:  {summary.final_clock_bias_error_m * 1e2:8.3f} cm")
    log(f"Final clock drift error: {summary.final_clock_drift_error_mps * 1e3:8.3f} mm/s")
    log(f"RMS position error:      {summary.rms_position_error_m:8.3f} m")
    log(f"RMS velocity error:      {summary.rms_velocity_error_mps * 1e3:8.3f} mm/s")
    if app_cfg.estimate_srp_coefficient:
        srp_est = app_cfg.srp_coeff_truth_m2_kg + summary.final_srp_coeff_error_m2_kg
        log(
            f"Final SRP coeff est:     {srp_est * 1e3:8.3f} e-3 m^2/kg "
            f"(truth {app_cfg.srp_coeff_truth_m2_kg * 1e3:.3f}e-3, "
            f"err {summary.final_srp_coeff_error_m2_kg * 1e3:+.3f}e-3)"
        )


def main():
    parser = argparse.ArgumentParser(description="Run the agent-based ex6 lunar GNSS ODTS EKF.")
    parser.add_argument(
        "--progress-epochs",
        type=int,
        default=None,
        help="Override run_progress_interval_epochs (0 = C++ default cadence).",
    )
    parser.add_argument(
        "--require-cache",
        action="store_true",
        help="Require the Stage 1/2 precompute cache to exist (else the EKF "
        "computes links in-memory, single-threaded).",
    )
    args = parser.parse_args()

    # ex6_gnss_odts_config sets LUPNT_DATA_PATH + sys.path before importing pylupnt.
    from ex6_gnss_odts_config import DATA_DIR, OUTPUT_DIR, load_scenario, describe
    import pylupnt as pnt

    scen = load_scenario()
    app_block = scen["agents"]["receiver"]["application"]
    if args.progress_epochs is not None:
        if args.progress_epochs < 0:
            parser.error("--progress-epochs must be >= 0")
        app_block["simulation"]["run_progress_interval_epochs"] = args.progress_epochs

    links = Path(app_block["pipeline"]["links_file"])
    delays = Path(app_block["pipeline"]["delays_file"])
    if args.require_cache:
        log("Checking precomputed inputs...")
        try:
            _require_input_file(links, "Stage 1 links CSV")
            if app_block["plasma"].get("simulate_truth"):
                _require_input_file(delays, "Stage 2 plasma-delay CSV")
        except FileNotFoundError as exc:
            print(str(exc), file=sys.stderr, flush=True)
            return 2

    log(f"Data:   {DATA_DIR}")
    log(f"Output: {OUTPUT_DIR}")
    log()

    sim = pnt.Simulation(scen)
    log("Starting sim.run() (agent-based LunarGnssOdtsApp)...")
    t0 = time.time()
    sim.run()
    log(f"sim.run() took {time.time() - t0:.1f} s")
    log()

    app = sim.get_agent("receiver").get_application()
    summaries = app.get_summaries()
    app_cfg = app.get_config()
    if not summaries:
        print("sim.run() completed but returned no summaries.", file=sys.stderr, flush=True)
        return 1
    _print_summary(summaries[0], app_cfg)

    # Per-seed summary export for programmatic access (the full time series live in the CSVs).
    npz_path = OUTPUT_DIR / "ex6_results.npz"
    np.savez(
        npz_path,
        final_position_error_m=np.array([s.final_position_error_m for s in summaries]),
        final_velocity_error_mps=np.array([s.final_velocity_error_mps for s in summaries]),
        final_clock_bias_error_m=np.array([s.final_clock_bias_error_m for s in summaries]),
        final_clock_drift_error_mps=np.array([s.final_clock_drift_error_mps for s in summaries]),
        final_srp_coeff_error_m2_kg=np.array([s.final_srp_coeff_error_m2_kg for s in summaries]),
        rms_position_error_m=np.array([s.rms_position_error_m for s in summaries]),
        rms_velocity_error_mps=np.array([s.rms_velocity_error_mps for s in summaries]),
        num_epochs=np.array([s.num_epochs for s in summaries]),
    )
    log()
    traj_path = OUTPUT_DIR / "trajectory_mc0.csv"
    if traj_path.exists():
        traj = pd.read_csv(traj_path)
        log(f"Wrote trajectory: {traj_path} ({len(traj)} rows)")
    log(f"Wrote summaries:  {npz_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
