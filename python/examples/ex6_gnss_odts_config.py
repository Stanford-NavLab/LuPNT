"""Shared configuration for the ex6 lunar GNSS ODTS example.

Imported by BOTH the notebook (``ex6_gnss_odts.ipynb``) and the standalone precompute
script (``ex6_precompute.py``) so they share the same baseline ``LunarGnssODTSConfig``. The C++
precompute cache is keyed only on settings that affect generated link geometry and truth
measurement metadata, so filter-only tuning and measurement-combination switches can reuse the
same cached measurements.
"""

import os
import sys
from pathlib import Path

# --- Use THIS repo's pylupnt/data, not a stale copy that a global PYTHONPATH (e.g. an
# `export PYTHONPATH=.../other/python` in ~/.zshrc) or a reverted Jupyter kernelspec might
# place first. This module knows the repo from its own location. ---
_REPO = Path(__file__).resolve().parents[2]  # python/examples/ -> repo root
_PYTHON = str((_REPO / "python").resolve())
if _PYTHON in sys.path:
    sys.path.remove(_PYTHON)
sys.path.insert(0, _PYTHON)

DATA_DIR = (_REPO / "data" / "LuPNT_data").resolve()
# Force LUPNT_DATA_PATH to this repo's data BEFORE importing pylupnt. pylupnt runs a
# download-on-import (pylupnt/core/download_data.py) that dumps a fresh LuPNT_data/ into the
# current working directory unless LUPNT_DATA_PATH already points at a dir containing
# `ephemeris/`. Setting it here (not setdefault -- override any stale value) keeps the data
# at the repo root and stops that stray download.
if (DATA_DIR / "ephemeris").is_dir():
    os.environ["LUPNT_DATA_PATH"] = str(DATA_DIR)

import numpy as np  # noqa: E402
import pylupnt as pnt  # noqa: E402

# Cache/output dir, anchored at the repo root (output/python_examples/) so the notebook
# (cwd = examples/) and the script (any cwd) always agree on the location.
OUTPUT_DIR = _REPO / "output" / "python_examples" / "ex6_gnss_odts_data"

# Ray-trace the (slow) plasmasphere delay once every 2 minutes and interpolate to 1 Hz.
PLASMA_DELAY_DT_S = 120.0

# Print Stage 2 ray-trace progress every N completed rays.
PLASMA_DELAY_PROGRESS_RAYS = 10


def build_config():
    """Build the ex6 ``LunarGnssODTSConfig`` (identical for the notebook and the script)."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    cfg = pnt.LunarGnssODTSConfig()

    # --- Timing: one full ELFO orbital period, processed at 1 Hz ---
    cfg.seed = 42
    cfg.monte_carlo_runs = 1
    cfg.start_epoch_utc = (
        "2026-01-01T02:00:00"  # COD MGEX SP3 (Galileo); offset fits its 1-day span
    )
    cfg.dt_s = 1.0
    cfg.ephemeris_dt_s = 30.0
    cfg.receiver_app.rate_hz = 1.0
    cfg.output_dir = str(OUTPUT_DIR)
    cfg.links_file = str(OUTPUT_DIR / "precomputed_links.csv")
    cfg.delays_file = str(OUTPUT_DIR / "precomputed_delays.csv")

    # --- Truth receiver orbit: ELFO (a=6541.4 km, e=0.6) ---
    cfg.receiver_a_m = 6541.4e3
    cfg.receiver_ecc = 0.6
    cfg.receiver_inc_rad = 65.5 * pnt.RAD
    cfg.receiver_raan_rad = 60.0 * pnt.RAD
    cfg.receiver_argp_rad = 90.0 * pnt.RAD
    cfg.receiver_mean_anomaly_rad = 0.0
    cfg.clock_bias_s = 2.0e-6
    cfg.clock_drift_sps = 1.0e-10
    # Receiver clock model (truth + filter). Options: OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC
    # (noise coefficients per model in ClockDynamics::GetClockValues).
    cfg.clock_model = "OCXO"
    # Optionally generate the *truth* clock with a 3-state [bias, drift, drift-rate] model
    # while the filter keeps its 2-state model (unmodeled clock aging). Off by default.
    cfg.use_three_state_clock_truth = False
    cfg.clock_drift_rate_sps2 = 0.0  # initial truth drift-rate [s/s^2] (3-state truth only)
    cfg.duration_s = float(2.0 * np.pi * np.sqrt(cfg.receiver_a_m**3 / pnt.GM_MOON))

    # --- GPS + Galileo constellation; auto-select the covering SP3 for the epoch ---
    sp3_dir = (DATA_DIR / "ephemeris" / "gnsslibpy" / "sp3").resolve()
    cfg.constellation.sp3_directory = str(sp3_dir)
    cfg.constellation.antex_file = str((DATA_DIR / "gnss" / "igs20.atx").resolve())
    cfg.constellation.auto_select_sp3 = True
    cfg.constellation.use_all_gps = True
    cfg.constellation.include_galileo = True

    # --- Broadcast (BRDC) transmitter ephemeris for the receiver model ---
    # Truth keeps the precise SP3 transmitter states; the filter (receiver) measurement model
    # is fed the broadcast position/clock, evaluated live from the RINEX-nav parameters at each
    # transmit epoch. The injected broadcast-minus-precise error is debiased first: the per-
    # constellation systematic clock offset and (for QZSS) the per-satellite radial orbit offset
    # are removed (Montenbruck & Steigenberger, J. Navigation, 2018). BRDC files must cover the
    # epoch (download alongside the SP3, e.g. via pylupnt.interfaces.gnss_file_loader.BRDCLoader).
    cfg.constellation.brdc_directory = str((DATA_DIR / "ephemeris" / "gnsslibpy" / "brdc").resolve())
    cfg.constellation.use_broadcast_ephemeris = True
    cfg.constellation.debias_broadcast_clock = True
    cfg.constellation.debias_qzss_radial = True

    # --- Sidelobe link budget ---
    cfg.design.receiver_params.Bp = 1.0
    cfg.design.receiver_params.T = 0.02
    cfg.design.receiver_params.b = 2.0
    cfg.design.receiver_params.Bn = 0.7
    cfg.design.receiver_params.Bf = 0.2
    cfg.design.receiver_params.D = 0.1
    cfg.design.setup_transmitters = True
    cfg.design.receiver_antenna_name = "moongpsr"
    cfg.design.apply_cn0_threshold = True
    cfg.design.cn0_acquisition_threshold_dbhz = 22.0
    cfg.design.cn0_tracking_threshold_dbhz = 20.0
    cfg.design.use_cn0_measurement_sigmas = True

    # --- Force model: 20x20 Moon gravity truth / 18x18 filter, Earth+Sun, relativity ---
    cfg.moon_gravity_degree_truth, cfg.moon_gravity_order_truth = 20, 20
    cfg.moon_gravity_degree_filter, cfg.moon_gravity_order_filter = 18, 18
    cfg.include_earth = True
    cfg.include_sun = True
    cfg.use_relativity = True

    # --- Solar radiation pressure: perturb truth and estimate the coefficient ---
    cfg.use_srp_truth = True
    cfg.srp_coeff_truth_m2_kg = 2.0e-3
    cfg.use_srp_filter = True
    cfg.estimate_srp_coefficient = True
    cfg.initial_srp_coeff_m2_kg = 1.0e-3
    cfg.initial_srp_coeff_sigma_m2_kg = 1.0e-3
    # Small random-walk floor to keep the SRP-coefficient variance from collapsing (numerical).
    cfg.process_srp_coeff_sigma_m2_kg_sqrt_s = 1.0e-9

    # --- Measurements: ionosphere-free (L1+L5) pseudorange + TDCP (L1), C/N0 noise ---
    cfg.use_pseudorange = True
    cfg.use_doppler = True
    cfg.use_tdcp = True
    cfg.use_ionosphere_free = True
    cfg.filter_pseudorange_noise_inflation_m = 10.0
    cfg.filter_tdcp_noise_inflation_m = 0.05
    cfg.pseudorange_min_tangent_altitude_m = 1000.0
    cfg.tdcp_min_tangent_altitude_m = 4000.0

    # --- Plasmaspheric delay: GCPM/IRI2007 ray trace, truth-only ---
    cfg.plasma.simulate_truth = True
    cfg.plasma.model_in_filter = False
    cfg.plasma.raytrace_kp = 3.0  # geomagnetic Kp index
    cfg.plasma.raytrace_rz12 = 50.0  # IRI R12 sunspot index (-1 = historical/projected)

    # --- Filter tuning ---
    cfg.initial_position_sigma_m = 100.0
    cfg.initial_velocity_sigma_mps = 0.1
    cfg.initial_clock_bias_sigma_s = 1.0e-6
    cfg.initial_clock_drift_sigma_sps = 1.0e-11
    cfg.process_accel_sigma_mps2 = 1.0e-9
    cfg.integration_step_s = 20.0

    # --- Diagnostics / performance ---
    # Minimum wall-clock seconds between Stage 1 precompute progress prints. Raise it to
    # reduce log spam on the long full-orbit run; lower it for more frequent updates.
    cfg.precompute_progress_interval_s = 5.0
    # Epochs between sim.run() EKF progress prints (position/clock error, tracked sats,
    # measurements). 0 keeps the C++ default of about 50 updates across the run.
    cfg.run_progress_interval_epochs = 1000
    # Print the STM and first few measurement-Jacobian rows for early filter epochs.
    cfg.debug_print_matrix_epochs = 3
    cfg.debug_print_matrix_max_rows = 8
    # Threads for the Stage 1 constellation loop (GPS L1/L5 + Galileo E1/E5a run in parallel).
    # 0 = all cores; the effective count is capped to the number of constellations.
    cfg.precompute_num_threads = 0

    return cfg


def describe(cfg):
    """One-line-per-topic summary of the config (used by the notebook and script)."""
    lines = [
        f"Configured {cfg.duration_s / 3600.0:.2f}-h run (~1 ELFO orbital period, "
        f"{cfg.duration_s:.0f} s) at {cfg.receiver_app.rate_hz:.0f} Hz "
        f"starting {cfg.start_epoch_utc}",
        f"  TDCP {'on' if cfg.use_tdcp else 'off'}, "
        f"SRP-coeff estimation {'on' if cfg.estimate_srp_coefficient else 'off'} "
        f"(truth Cr*A/m = {cfg.srp_coeff_truth_m2_kg * 1e3:.2f}e-3 m^2/kg)",
        f"  Pseudorange: {'ionosphere-free (L1+L5)' if cfg.use_ionosphere_free else 'L1'}, "
        f"C/N0 noise + filter inflation (PR {cfg.filter_pseudorange_noise_inflation_m:.0f} m, "
        f"TDCP {cfg.filter_tdcp_noise_inflation_m * 1e3:.0f} mm)",
        f"  Tangent-altitude cutoff: PR >= {cfg.pseudorange_min_tangent_altitude_m:.0f} m, "
        f"TDCP >= {cfg.tdcp_min_tangent_altitude_m:.0f} m",
        f"  Constellation: GPS{' + Galileo' if cfg.constellation.include_galileo else ''}, "
        f"receiver antenna {cfg.design.receiver_antenna_name or 'omni'}, "
        f"C/N0 gate {cfg.design.cn0_acquisition_threshold_dbhz:.0f}/"
        f"{cfg.design.cn0_tracking_threshold_dbhz:.0f} dB-Hz (acquisition/tracking)",
        f"  Plasma delay ray-traced every {PLASMA_DELAY_DT_S:.0f} s, interpolated to "
        f"{cfg.receiver_app.rate_hz:.0f} Hz; progress every {PLASMA_DELAY_PROGRESS_RAYS} rays",
        f"  EKF run progress every {cfg.run_progress_interval_epochs} epochs "
        f"({cfg.run_progress_interval_epochs / cfg.receiver_app.rate_hz:.0f} s)",
        f"Data directory: {OUTPUT_DIR}",
    ]
    return "\n".join(lines)
