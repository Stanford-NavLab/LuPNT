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

import yaml  # noqa: E402
import pylupnt as pnt  # noqa: E402

# Cache/output dir, anchored at the repo root (output/python_examples/) so the notebook
# (cwd = examples/) and the script (any cwd) always agree on the location.
OUTPUT_DIR = _REPO / "output" / "python_examples" / "ex6_gnss_odts_data"

# Agent-based scenario config (world: + a physical `receiver` Spacecraft hosting a LunarGnssOdtsApp). The
# EKF run is driven from this via ``pnt.Simulation``; its scalar values mirror build_config()
# exactly so the C++ link-precompute cache (keyed on measurement-generation inputs) is shared.
CONFIG_YAML = _REPO / "configs" / "lunar_gnss_odts.yaml"

# Ray-trace the (slow) plasmasphere delay once every 2 minutes and interpolate to 1 Hz.
PLASMA_DELAY_DT_S = 120.0

# Print Stage 2 ray-trace progress every N completed rays.
PLASMA_DELAY_PROGRESS_RAYS = 10


def _truth_from_yaml():
    """Read the receiver's TRUTH (orbit + clock + force model) and the run timing from
    ``configs/lunar_gnss_odts.yaml`` -- the single source of truth.

    In the agent-based scenario this truth is declared once, on the physical ``receiver``
    Spacecraft (its ``initial_state:``/``dynamics:`` blocks) and the app ``simulation:`` block.
    Deriving build_config()'s struct fields from the same YAML (instead of a hand-kept copy)
    means the standalone precompute and the agent EKF run can never drift out of the shared
    link-cache fingerprint.
    """
    scen = yaml.safe_load(open(CONFIG_YAML))
    rx = scen["agents"]["receiver"]
    init, dyn = rx["initial_state"], rx["dynamics"]
    sim = rx["application"]["simulation"]
    # bodies: a list of single-key maps, e.g. [{MOON: {n: 20, m: 20}}, {EARTH: {}}, {SUN: {}}]
    bodies = {next(iter(b)): (b[next(iter(b))] or {}) for b in dyn.get("bodies", [])}
    moon = bodies.get("MOON", {})
    mass, area, cr = float(dyn.get("mass", 1.0)), float(dyn.get("area", 0.0)), float(dyn.get("CR", 0.0))
    return dict(
        start_epoch_utc=scen.get("epoch", sim.get("start_epoch_utc")),
        seed=int(sim.get("seed", 42)),
        dt_s=float(sim["dt_s"]),
        ephemeris_dt_s=float(sim["ephemeris_dt_s"]),
        receiver_rate_hz=float(rx["application"]["receiver_app"]["rate_hz"]),
        duration_s=float(sim["duration_s"]),
        receiver_a_m=float(init["a"]),
        receiver_ecc=float(init["e"]),
        receiver_inc_rad=float(init["i"]) * pnt.RAD,
        receiver_raan_rad=float(init["Omega"]) * pnt.RAD,
        receiver_argp_rad=float(init["omega"]) * pnt.RAD,
        receiver_mean_anomaly_rad=float(init["M"]) * pnt.RAD,
        clock_bias_s=float(init.get("clock_bias_s", 0.0)),
        clock_drift_sps=float(init.get("clock_drift_sps", 0.0)),
        moon_gravity_degree_truth=int(moon.get("n", 0)),
        moon_gravity_order_truth=int(moon.get("m", 0)),
        include_earth="EARTH" in bodies,
        include_sun="SUN" in bodies,
        use_relativity=bool(dyn.get("use_relativity", False)),
        use_srp_truth=area > 0.0,
        srp_coeff_truth_m2_kg=(cr * area / mass) if mass else 0.0,
    )


def build_config():
    """Build the ex6 ``LunarGnssODTSConfig`` (identical for the notebook and the script).

    The truth orbit/clock/force-model + run timing are read from ``configs/lunar_gnss_odts.yaml``
    (see ``_truth_from_yaml``) so this struct config and the agent-run config share one source and
    stay fingerprint-compatible; only the constellation, sidelobe, and filter-only tunings are set
    here.
    """
    OUTPUT_DIR.mkdir(exist_ok=True)
    cfg = pnt.LunarGnssODTSConfig()
    _T = _truth_from_yaml()

    # --- Timing + truth receiver orbit/clock: read from the YAML (single source) ---
    cfg.monte_carlo_runs = 1
    cfg.start_epoch_utc = _T["start_epoch_utc"]  # COD MGEX SP3 (Galileo); offset fits its 1-day span
    cfg.seed = _T["seed"]
    cfg.dt_s = _T["dt_s"]
    cfg.ephemeris_dt_s = _T["ephemeris_dt_s"]
    cfg.receiver_app.rate_hz = _T["receiver_rate_hz"]
    cfg.duration_s = _T["duration_s"]  # one full ELFO orbital period
    cfg.output_dir = str(OUTPUT_DIR)
    cfg.links_file = str(OUTPUT_DIR / "precomputed_links.csv")
    cfg.delays_file = str(OUTPUT_DIR / "precomputed_delays.csv")

    # Truth receiver ELFO orbit + clock (from the `receiver` Spacecraft's initial_state).
    cfg.receiver_a_m = _T["receiver_a_m"]
    cfg.receiver_ecc = _T["receiver_ecc"]
    cfg.receiver_inc_rad = _T["receiver_inc_rad"]
    cfg.receiver_raan_rad = _T["receiver_raan_rad"]
    cfg.receiver_argp_rad = _T["receiver_argp_rad"]
    cfg.receiver_mean_anomaly_rad = _T["receiver_mean_anomaly_rad"]
    cfg.clock_bias_s = _T["clock_bias_s"]
    cfg.clock_drift_sps = _T["clock_drift_sps"]
    # Receiver clock model (truth + filter). Options: OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC
    # (noise coefficients per model in ClockDynamics::GetClockValues).
    cfg.clock_model = "OCXO"
    # Optionally generate the *truth* clock with a 3-state [bias, drift, drift-rate] model
    # while the filter keeps its 2-state model (unmodeled clock aging). Off by default.
    cfg.use_three_state_clock_truth = False
    cfg.clock_drift_rate_sps2 = 0.0  # initial truth drift-rate [s/s^2] (3-state truth only)

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
    # epoch (download alongside the SP3, e.g. via pnt.RinexNavLoader.download_file_for_epoch).
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
    # Transmitter yaw-steering model for the C/N0 antenna-gain geometry. False (default) uses the
    # nominal Sun-pointing frame; set True to use the block-specific *dedicated* eclipse
    # yaw-steering law (GPS/Galileo) -- it differs from nominal only near orbit noon/midnight at
    # low Sun-beta, where sidelobe reception is most affected. See ex3 for the attitude laws.
    cfg.design.tx_yaw_dedicated = False

    # --- Force model: Moon gravity truth (from the receiver Spacecraft's dynamics), 18x18
    #     filter, Earth+Sun, relativity ---
    cfg.moon_gravity_degree_truth = _T["moon_gravity_degree_truth"]
    cfg.moon_gravity_order_truth = _T["moon_gravity_order_truth"]
    cfg.moon_gravity_degree_filter, cfg.moon_gravity_order_filter = 18, 18
    cfg.include_earth = _T["include_earth"]
    cfg.include_sun = _T["include_sun"]
    cfg.use_relativity = _T["use_relativity"]

    # --- Solar radiation pressure: perturb truth (Cr*A/m from the Spacecraft dynamics) and
    #     estimate the coefficient ---
    cfg.use_srp_truth = _T["use_srp_truth"]
    cfg.srp_coeff_truth_m2_kg = _T["srp_coeff_truth_m2_kg"]
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


def load_scenario():
    """Load ``configs/lunar_gnss_odts.yaml`` and patch its data/output paths to THIS repo.

    Returns the full scenario dict (top-level ``world:`` + a physical ``receiver`` Spacecraft agent hosting a
    ``LunarGnssOdtsApp``) ready to pass to ``pnt.Simulation``. The application-block scalar
    values match ``build_config()`` field-for-field, so an EKF run driven from this dict reuses
    the Stage 1/2 precompute cache generated by ``ex6_precompute.py`` (build_config path).
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_YAML) as f:
        scen = yaml.safe_load(f)
    app = scen["agents"]["receiver"]["application"]
    app["simulation"]["output_dir"] = str(OUTPUT_DIR)
    app["pipeline"]["links_file"] = str(OUTPUT_DIR / "precomputed_links.csv")
    app["pipeline"]["delays_file"] = str(OUTPUT_DIR / "precomputed_delays.csv")
    c = app["constellation"]
    c["sp3_directory"] = str((DATA_DIR / "ephemeris" / "gnsslibpy" / "sp3").resolve())
    c["antex_file"] = str((DATA_DIR / "gnss" / "igs20.atx").resolve())
    c["brdc_directory"] = str((DATA_DIR / "ephemeris" / "gnsslibpy" / "brdc").resolve())

    # --- Link-cache fingerprint compatibility -------------------------------------------------
    # The Stage 1/2 precompute (ex6_precompute.py) keys its link cache on a fingerprint of the
    # TRUTH link geometry -- the receiver orbit + truth force model -- computed from build_config().
    # In the agent-based scenario that truth is owned by the physical `receiver` Spacecraft (its
    # `dynamics:`/`initial_state:` blocks), so the *application* block deliberately does not declare
    # it. But the C++ fingerprint still reads a `truth:`/`dynamics:` block, so we mirror
    # build_config()'s exact values here (identical doubles, same source) to reconstruct the cache
    # key so the EKF run reuses the precompute cache. These do NOT drive the run -- the truth comes
    # from the Spacecraft agent; they only reproduce the fingerprint.
    ref = build_config()
    t = app.setdefault("truth", {})
    t["receiver_a_m"] = ref.receiver_a_m
    t["receiver_ecc"] = ref.receiver_ecc
    t["receiver_inc_rad"] = ref.receiver_inc_rad
    t["receiver_raan_rad"] = ref.receiver_raan_rad
    t["receiver_argp_rad"] = ref.receiver_argp_rad
    t["receiver_mean_anomaly_rad"] = ref.receiver_mean_anomaly_rad
    t["clock_bias_s"] = ref.clock_bias_s
    t["clock_drift_sps"] = ref.clock_drift_sps
    d = app.setdefault("dynamics", {})
    d["moon_gravity_degree_truth"] = ref.moon_gravity_degree_truth
    d["moon_gravity_order_truth"] = ref.moon_gravity_order_truth
    d["include_earth"] = ref.include_earth
    d["include_sun"] = ref.include_sun
    d["use_relativity"] = ref.use_relativity
    d["use_srp_truth"] = ref.use_srp_truth
    d["srp_coeff_truth_m2_kg"] = ref.srp_coeff_truth_m2_kg
    return scen


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
