"""Run the ex6 GNSS ODTS Kalman-filter stage from a terminal.

This is the standalone equivalent of the notebook cell that calls ``sim.run()`` in
``ex6_gnss_odts.ipynb``. Run the slower precompute stages first:

    pixi run python python/examples/ex6_precompute.py

Then run this script:

    pixi run python python/examples/ex6_run_kalman.py

All Python prints use ``flush=True`` so progress is visible promptly when the script is run
from a terminal or launched from another process.
"""

import argparse
import sys
import time
from pathlib import Path

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


def _require_valid_link_cache(cfg, pnt):
    if pnt.lunar_gnss_odts_link_cache_valid(cfg):
        return
    raise RuntimeError(
        "The Stage 1 link cache exists, but it is stale or was written with an older CSV "
        "schema. Regenerate Stage 1/2 before running the EKF:\n"
        "  pixi run python python/examples/ex6_precompute.py --force"
    )


def _remove_previous_outputs(cfg):
    output_dir = Path(cfg.output_dir)
    removed = []
    for path in sorted(output_dir.glob("trajectory_mc*.csv")):
        path.unlink()
        removed.append(path)
    summary = output_dir / "summary.csv"
    if summary.exists():
        summary.unlink()
        removed.append(summary)
    return removed


def _progress_description(cfg):
    cadence = int(getattr(cfg, "run_progress_interval_epochs", 0))
    if cadence > 0:
        seconds = cadence / float(cfg.receiver_app.rate_hz)
        return f"every {cadence} epochs ({seconds:.0f} s)"
    return "C++ default cadence (~50 updates across the run)"


def _print_summary(summary, cfg):
    log(f"Final position error:    {summary.final_position_error_m:8.3f} m")
    log(f"Final velocity error:    {summary.final_velocity_error_mps * 1e3:8.3f} mm/s")
    log(f"Final clock bias error:  {summary.final_clock_bias_error_m * 1e2:8.3f} cm")
    log(f"Final clock drift error: {summary.final_clock_drift_error_mps * 1e3:8.3f} mm/s")
    log(f"RMS position error:      {summary.rms_position_error_m:8.3f} m")
    log(f"RMS velocity error:      {summary.rms_velocity_error_mps * 1e3:8.3f} mm/s")
    if cfg.estimate_srp_coefficient:
        srp_est = cfg.srp_coeff_truth_m2_kg + summary.final_srp_coeff_error_m2_kg
        log(
            f"Final SRP coeff est:     {srp_est * 1e3:8.3f} e-3 m^2/kg "
            f"(truth {cfg.srp_coeff_truth_m2_kg * 1e3:.3f}e-3, "
            f"err {summary.final_srp_coeff_error_m2_kg * 1e3:+.3f}e-3)"
        )


def _load_config_and_bindings():
    try:
        from ex6_gnss_odts_config import build_config, describe
        import pylupnt as pnt
    except ModuleNotFoundError as exc:
        if exc.name in {"pylupnt._pylupnt", "_pylupnt"}:
            raise SystemExit(
                "Could not import the compiled pylupnt extension.\n"
                "Rebuild the Python bindings for the active Pixi interpreter:\n"
                "  pixi run build-py"
            ) from exc
        raise
    return build_config, describe, pnt


def main():
    parser = argparse.ArgumentParser(
        description="Run Step 5 / EKF execution for python/examples/ex6_gnss_odts.ipynb."
    )
    parser.add_argument(
        "--progress-epochs",
        type=int,
        default=None,
        help=(
            "Override cfg.run_progress_interval_epochs for sim.run() progress prints. "
            "Use 0 for the C++ default cadence."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete previous trajectory_mc*.csv and summary.csv outputs before running.",
    )
    args = parser.parse_args()

    build_config, describe, pnt = _load_config_and_bindings()
    cfg = build_config()
    if args.progress_epochs is not None:
        if args.progress_epochs < 0:
            parser.error("--progress-epochs must be >= 0")
        cfg.run_progress_interval_epochs = args.progress_epochs

    log(describe(cfg))
    log()
    log("Checking precomputed inputs...")
    try:
        _require_input_file(cfg.links_file, "Stage 1 links CSV")
        _require_input_file(cfg.delays_file, "Stage 2 plasma-delay CSV")
        _require_valid_link_cache(cfg, pnt)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr, flush=True)
        return 2
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr, flush=True)
        return 2

    if args.force:
        removed = _remove_previous_outputs(cfg)
        if removed:
            log(f"Removed {len(removed)} previous run output file(s).")

    log(f"Links:  {cfg.links_file}")
    log(f"Delays: {cfg.delays_file}")
    log(f"Output: {cfg.output_dir}")
    log(f"EKF progress cadence: {_progress_description(cfg)}")
    log()

    sim = pnt.LunarGnssODTSSimulation(cfg)
    log("Starting sim.run()...")
    t0 = time.time()
    sim.run()
    elapsed_s = time.time() - t0
    log(f"sim.run() took {elapsed_s:.1f} s")
    log()

    summaries = sim.get_summaries()
    if not summaries:
        print("sim.run() completed but returned no summaries.", file=sys.stderr, flush=True)
        return 1
    _print_summary(summaries[0], cfg)

    traj_path = Path(cfg.output_dir) / "trajectory_mc0.csv"
    summary_path = Path(cfg.output_dir) / "summary.csv"
    if traj_path.exists():
        traj = pd.read_csv(traj_path)
        log()
        log(f"Wrote trajectory: {traj_path} ({len(traj)} rows)")
        if len(traj):
            final_t_hr = float(traj.iloc[-1]["t"]) / 3600.0
            log(f"Final trajectory time: {final_t_hr:.3f} h")
    if summary_path.exists():
        log(f"Wrote summary:    {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
