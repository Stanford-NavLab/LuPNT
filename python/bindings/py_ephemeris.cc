/**
 * @file py_ephemeris.cc
 * @brief Python bindings for `lupnt::LansEphemeris`/`lupnt::LansAlmanac`
 *        (`lupnt/applications/ephemeris/lunanet_ephemeris.h`,
 * `lupnt/applications/ephemeris/lunanet_almanac.h`), the ephemeris/almanac datasize-accuracy study
 * config/result structs
 *        (`lupnt/simulations/ephemeris/ephemeris_simulation.h`), and the agent-based
 *        `lupnt::EphemerisApp` (`lupnt/applications/ephemeris/ephemeris_app.h`).
 *
 * Config/result/options structs are plain `double`-valued (no `lupnt::Real`), so
 * they bind directly via `def_readwrite`/`def_readonly` with pybind11's built-in
 * Eigen<->NumPy conversion, following the `IslOdtsConfig`/`IslOdtsResults` pattern
 * in `py_isl_odts.cc`.
 */
#include <lupnt/applications/ephemeris/ephemeris_app.h>
#include <lupnt/lupnt.h>

#include <memory>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitEphemeris(py::module& m) {
  // ---- EphemerisFitErrorStats -------------------------------------------------

  py::class_<EphemerisFitErrorStats>(
      m, "EphemerisFitErrorStats",
      "Position/velocity fit-error statistics in Radial/Transverse/Normal (RTN) components plus "
      "the combined 3-D norm.")
      .def(py::init<>(), "Default-construct (all errors zero).")
      .def_readonly("rms_pos_m", &EphemerisFitErrorStats::rms_pos_m,
                    "[R, T, N, 3D] position RMS error [m]")
      .def_readonly("rms_vel_mps", &EphemerisFitErrorStats::rms_vel_mps,
                    "[R, T, N, 3D] velocity RMS error [m/s]")
      .def_readonly("p95_pos_m", &EphemerisFitErrorStats::p95_pos_m,
                    "[R, T, N, 3D] 95th-percentile |position error| [m]")
      .def_readonly("p95_vel_mps", &EphemerisFitErrorStats::p95_vel_mps,
                    "[R, T, N, 3D] 95th-percentile |velocity error| [m/s]");

  // ---- LansEphemeris ------------------------------------------------------

  py::class_<EphemerisFitOptions>(m, "EphemerisFitOptions",
                                  "Configuration for LansEphemeris.fit/eval.")
      .def(py::init<>(), "Default-construct with default fit options.")
      .def_readwrite("poly_order", &EphemerisFitOptions::poly_order,
                     "Chebyshev polynomial order for the position residual (velocity from its "
                     "analytic derivative).")
      .def_readwrite("use_keplerian_baseline", &EphemerisFitOptions::use_keplerian_baseline,
                     "If true, subtract an osculating two-body Kepler baseline before the "
                     "polynomial fit.")
      .def_readwrite("gm", &EphemerisFitOptions::gm,
                     "Central-body gravitational parameter [m^3/s^2].")
      .def_readwrite("frame", &EphemerisFitOptions::frame,
                     "Frame the input states and fitted output are represented in (e.g. MOON_CI, "
                     "MOON_PA).")
      .def_readwrite("num_fourier_terms", &EphemerisFitOptions::num_fourier_terms,
                     "Fourier terms (harmonics of the argument of latitude) added to the residual "
                     "model; 0 = pure Chebyshev.");

  py::class_<LansEphemeris>(
      m, "LansEphemeris",
      "Piecewise ephemeris model: an osculating two-body Kepler orbit plus a Chebyshev correction "
      "on the Cartesian position residual (GNSS broadcast-ephemeris style).")
      .def(py::init<>(), "Default-construct with default fit options.")
      .def(py::init<EphemerisFitOptions>(), py::arg("options"),
           "Construct with the given fit options.")
      .def("fit", &LansEphemeris::Fit, py::arg("t_s"), py::arg("rv"),
           "Fit the model to sampled states `rv` [N x 6] at epochs `t_s` [s]; returns the fitted "
           "parameter vector.")
      .def("eval", &LansEphemeris::Eval, py::arg("t_s"), py::arg("params"),
           "Evaluate the fitted ephemeris `params` at epochs `t_s` [s]; returns states [N x 6].")
      .def("eval_error", &LansEphemeris::EvalError, py::arg("t_s"), py::arg("rv_ref"),
           py::arg("params"),
           "Fit-error statistics of `params` against reference states `rv_ref` at epochs `t_s`.")
      .def("num_params", &LansEphemeris::NumParams,
           "Number of scalar parameters in the fitted vector.")
      .def("param_names", &LansEphemeris::ParamNames,
           "Human-readable name of each fitted parameter (same order as fit()).")
      .def("get_options", &LansEphemeris::GetOptions, py::return_value_policy::reference_internal,
           "The EphemerisFitOptions this model was constructed with.");

  // ---- LansAlmanac ------------------------------------------------------------------

  py::class_<AlmanacFitOptions>(m, "AlmanacFitOptions", "Configuration for LansAlmanac.fit/eval.")
      .def(py::init<>(), "Default-construct with default fit options.")
      .def_readwrite("poly_order", &AlmanacFitOptions::poly_order,
                     "Secular polynomial degree fit to each element's time history (default 1 = "
                     "linear).")
      .def_readwrite("num_fourier_terms", &AlmanacFitOptions::num_fourier_terms,
                     "Fourier harmonics added per element (element-specific base frequency).")
      .def_readwrite("gm", &AlmanacFitOptions::gm,
                     "Central-body gravitational parameter [m^3/s^2].")
      .def_readwrite("sidereal_period_s", &AlmanacFitOptions::sidereal_period_s,
                     "Central-body sidereal rotation period [s], setting the non-SMA Fourier base "
                     "frequency.")
      .def_readwrite("frame", &AlmanacFitOptions::frame,
                     "Frame the input states and fitted output are represented in (e.g. MOON_CI, "
                     "MOON_PA).");

  py::class_<LansAlmanac>(
      m, "LansAlmanac",
      "Coarse almanac-style orbit model (low-order polynomial + Fourier per osculating element); "
      "GNSS-almanac style, long validity / small size.")
      .def(py::init<>(), "Default-construct with default fit options.")
      .def(py::init<AlmanacFitOptions>(), py::arg("options"),
           "Construct with the given fit options.")
      .def("fit", &LansAlmanac::Fit, py::arg("t_s"), py::arg("rv"),
           "Fit the model to sampled states `rv` [N x 6] at epochs `t_s` [s]; returns the fitted "
           "parameter vector.")
      .def("eval", &LansAlmanac::Eval, py::arg("t_s"), py::arg("params"),
           "Evaluate the fitted almanac `params` at epochs `t_s` [s]; returns states [N x 6].")
      .def("eval_error", &LansAlmanac::EvalError, py::arg("t_s"), py::arg("rv_ref"),
           py::arg("params"),
           "Fit-error statistics of `params` against reference states `rv_ref` at epochs `t_s`.")
      .def("num_params", &LansAlmanac::NumParams,
           "Number of scalar parameters in the fitted vector.")
      .def("param_names", &LansAlmanac::ParamNames,
           "Human-readable name of each fitted parameter (same order as fit()).")
      .def("get_options", &LansAlmanac::GetOptions, py::return_value_policy::reference_internal,
           "The AlmanacFitOptions this model was constructed with.");

  // ---- Ephemeris study config/results (EphemerisApp) -----------------------------

  py::class_<EphemerisOrbitConfig>(
      m, "EphemerisOrbitConfig",
      "Initial classical orbital elements for the study's truth trajectory.")
      .def(py::init<>(), "Default-construct with default orbital elements.")
      .def_readwrite("a_m", &EphemerisOrbitConfig::a_m, "Semi-major axis [m].")
      .def_readwrite("ecc", &EphemerisOrbitConfig::ecc, "Eccentricity.")
      .def_readwrite("inc_rad", &EphemerisOrbitConfig::inc_rad, "Inclination [rad].")
      .def_readwrite("raan_rad", &EphemerisOrbitConfig::raan_rad,
                     "Right ascension of the ascending node [rad].")
      .def_readwrite("argp_rad", &EphemerisOrbitConfig::argp_rad, "Argument of periapsis [rad].")
      .def_readwrite("m0_rad", &EphemerisOrbitConfig::m0_rad, "Initial mean anomaly [rad].")
      .def_readwrite("coe_frame", &EphemerisOrbitConfig::coe_frame,
                     "Frame the elements are defined in (e.g. MOON_OP).");

  py::class_<EphemerisSimulationConfig>(
      m, "EphemerisSimulationConfig",
      "Configuration for the ephemeris/almanac datasize-accuracy study.")
      .def(py::init<>(), "Default-construct with default study settings.")
      .def_readwrite("start_epoch_utc", &EphemerisSimulationConfig::start_epoch_utc,
                     "Start epoch (UTC ISO-8601 string).")
      .def_readwrite("orbit", &EphemerisSimulationConfig::orbit,
                     "Initial orbital elements of the truth trajectory.")
      .def_readwrite("propagate_frame", &EphemerisSimulationConfig::propagate_frame,
                     "Frame the truth trajectory is propagated and models are fit in (e.g. "
                     "MOON_CI).")
      .def_readwrite("duration_days", &EphemerisSimulationConfig::duration_days,
                     "Truth-trajectory duration [days].")
      .def_readwrite("sample_dt_s", &EphemerisSimulationConfig::sample_dt_s,
                     "Truth-trajectory sampling interval [s].")
      .def_readwrite("moon_gravity_degree", &EphemerisSimulationConfig::moon_gravity_degree,
                     "Moon gravity-field degree for the truth force model.")
      .def_readwrite("moon_gravity_order", &EphemerisSimulationConfig::moon_gravity_order,
                     "Moon gravity-field order for the truth force model.")
      .def_readwrite("include_earth", &EphemerisSimulationConfig::include_earth,
                     "Include Earth third-body gravity in the truth force model.")
      .def_readwrite("include_sun", &EphemerisSimulationConfig::include_sun,
                     "Include Sun third-body gravity in the truth force model.")
      .def_readwrite("use_relativity", &EphemerisSimulationConfig::use_relativity,
                     "Include relativistic corrections in the truth force model.")
      .def_readwrite("integration_step_s", &EphemerisSimulationConfig::integration_step_s,
                     "Truth-trajectory integration step [s].")
      .def_readwrite("fit_window_minutes", &EphemerisSimulationConfig::fit_window_minutes,
                     "Fitting-window lengths swept by the study [min].")
      .def_readwrite("num_windows", &EphemerisSimulationConfig::num_windows,
                     "Number of fitting windows sampled per fit_window_minutes entry.")
      .def_readwrite("cartesian_poly_order", &EphemerisSimulationConfig::cartesian_poly_order,
                     "Chebyshev polynomial order for the LansEphemeris model.")
      .def_readwrite("cartesian_use_keplerian_baseline",
                     &EphemerisSimulationConfig::cartesian_use_keplerian_baseline,
                     "Use a two-body Kepler baseline for the LansEphemeris model.")
      .def_readwrite("cartesian_num_fourier_terms",
                     &EphemerisSimulationConfig::cartesian_num_fourier_terms,
                     "Fourier terms for the LansEphemeris model.")
      .def_readwrite("almanac_poly_order", &EphemerisSimulationConfig::almanac_poly_order,
                     "Secular polynomial degree for the LansAlmanac model.")
      .def_readwrite("almanac_num_fourier_terms",
                     &EphemerisSimulationConfig::almanac_num_fourier_terms,
                     "Fourier harmonics per element for the LansAlmanac model.")
      .def_readwrite("output_frame", &EphemerisSimulationConfig::output_frame,
                     "Frame the models are fit and evaluated in (defaults to propagate_frame; "
                     "e.g. MOON_PA).")
      .def_readwrite("datasize_precision_m", &EphemerisSimulationConfig::datasize_precision_m,
                     "Required position accuracy used to size each parameter's broadcast "
                     "resolution [m].");

  py::class_<EphemerisWindowResult>(
      m, "EphemerisWindowResult",
      "Fit-accuracy / broadcast-datasize summary for one fit_window_minutes entry.")
      .def(py::init<>(), "Default-construct (all fields zero).")
      .def_readonly("fit_window_min", &EphemerisWindowResult::fit_window_min,
                    "Fitting-window length [min].")
      .def_readonly("num_params", &EphemerisWindowResult::num_params, "Number of model parameters.")
      .def_readonly("total_bits", &EphemerisWindowResult::total_bits,
                    "Total broadcast message size [bits].")
      .def_readonly("pos_rms_m", &EphemerisWindowResult::pos_rms_m, "Position RMS error [m].")
      .def_readonly("vel_rms_mps", &EphemerisWindowResult::vel_rms_mps, "Velocity RMS error [m/s].")
      .def_readonly("pos_p95_m", &EphemerisWindowResult::pos_p95_m,
                    "95th-percentile position error [m].")
      .def_readonly("vel_p95_mps", &EphemerisWindowResult::vel_p95_mps,
                    "95th-percentile velocity error [m/s].");

  py::class_<EphemerisResults>(
      m, "EphemerisResults",
      "Full result payload of an ephemeris datasize/accuracy sweep (SI units, output_frame).")
      .def(py::init<>(), "Default-construct (empty results).")
      .def_readonly("t_truth_s", &EphemerisResults::t_truth_s,
                    "Elapsed time since start_epoch_utc [s], size [N]")
      .def_readonly("rv_truth", &EphemerisResults::rv_truth,
                    "Truth Cartesian states [N x 6] in output_frame")
      .def_readonly("cartesian_results", &EphemerisResults::cartesian_results,
                    "List of EphemerisWindowResult, one per fit_window_minutes value "
                    "(LansEphemeris)")
      .def_readonly(
          "almanac_results", &EphemerisResults::almanac_results,
          "List of EphemerisWindowResult, one per fit_window_minutes value (LansAlmanac)");

  // ---- EphemerisApp (agent-based coordinator) --------------------------------------
  // Hosted on an `EphemerisManager` agent. Retrieve it from a `pnt.Simulation` via
  // `sim.get_agent("EphemerisManager").get_application()` (downcasts here), then read its
  // `EphemerisResults`.
  py::class_<EphemerisApp, Application, std::shared_ptr<EphemerisApp>>(
      m, "EphemerisApp",
      "Coordinator application for the ephemeris/almanac datasize-accuracy study (Example 9), "
      "hosted on an EphemerisManager agent.")
      .def(py::init<EphemerisSimulationConfig>(), py::arg("config"),
           "Construct from an EphemerisSimulationConfig.")
      .def("get_config", &EphemerisApp::GetConfig, py::return_value_policy::reference_internal,
           "The EphemerisSimulationConfig this app runs.")
      .def("get_results", &EphemerisApp::GetResults, py::return_value_policy::reference_internal,
           "Full EphemerisResults (truth trajectory + per-window Cartesian/LansAlmanac results)");

  // ---- EphemerisGenApp (LunaNet nav-message generation sub-app) ------------------

  py::class_<BroadcastMessage>(
      m, "BroadcastMessage",
      "One generated broadcast navigation message (ephemeris or almanac) valid over a single "
      "fitting window.")
      .def(py::init<>(), "Default-construct (empty message).")
      .def_readonly("t_generated_s", &BroadcastMessage::t_generated_s,
                    "Elapsed time the message was generated [s].")
      .def_readonly("t_start_s", &BroadcastMessage::t_start_s,
                    "Validity-window start [s, same origin].")
      .def_readonly("t_end_s", &BroadcastMessage::t_end_s, "Validity-window end [s, same origin].")
      .def_readonly("frame", &BroadcastMessage::frame,
                    "Frame of the fitted arc (= agent dynamics frame).")
      .def_readonly("params", &BroadcastMessage::params,
                    "Fitted parameter vector of the model (see the model's param_names()).");

  py::class_<EphemerisGenConfig>(m, "EphemerisGenConfig", "Configuration for EphemerisGenApp.")
      .def(py::init<>(), "Default-construct with default settings.")
      .def_readwrite("generate_ephemeris", &EphemerisGenConfig::generate_ephemeris,
                     "Generate a precise, short-validity broadcast ephemeris.")
      .def_readwrite("generate_almanac", &EphemerisGenConfig::generate_almanac,
                     "Generate a coarse, long-validity broadcast almanac.")
      .def_readwrite("ephemeris_options", &EphemerisGenConfig::ephemeris_options,
                     "Fit options for the broadcast ephemeris.")
      .def_readwrite("ephemeris_window_s", &EphemerisGenConfig::ephemeris_window_s,
                     "Ephemeris validity-window length [s].")
      .def_readwrite("ephemeris_refresh_s", &EphemerisGenConfig::ephemeris_refresh_s,
                     "Ephemeris refresh cadence [s] (<=0 = once per validity window).")
      .def_readwrite("almanac_options", &EphemerisGenConfig::almanac_options,
                     "Fit options for the broadcast almanac.")
      .def_readwrite("almanac_window_s", &EphemerisGenConfig::almanac_window_s,
                     "LansAlmanac validity-window length [s].")
      .def_readwrite("almanac_refresh_s", &EphemerisGenConfig::almanac_refresh_s,
                     "LansAlmanac refresh cadence [s] (<=0 = once per validity window).")
      .def_readwrite("ephemeris_fit_samples", &EphemerisGenConfig::ephemeris_fit_samples,
                     "Number of arc samples used to fit each ephemeris window.")
      .def_readwrite("almanac_fit_samples", &EphemerisGenConfig::almanac_fit_samples,
                     "Number of arc samples used to fit each almanac window.")
      .def_readwrite("output_frame", &EphemerisGenConfig::output_frame,
                     "Frame the broadcast messages are fit and evaluated in (e.g. MOON_PA).");

  py::class_<EphemerisGenApp>(
      m, "EphemerisGenApp",
      "LunaNet sub-app that generates broadcast ephemeris and almanac messages from the "
      "satellite's own predicted trajectory.")
      .def(py::init<>(), "Default-construct with default config.")
      .def(py::init<EphemerisGenConfig>(), py::arg("config"), "Construct with the given config.")
      .def("generate_ephemeris_from_arc", &EphemerisGenApp::GenerateEphemerisFromArc,
           py::arg("t_gen_s"), py::arg("t_s"), py::arg("rv"),
           py::return_value_policy::reference_internal,
           "Fit and store an ephemeris message from an explicit sampled arc (`t_s` [s], "
           "`rv` [N x 6]).")
      .def("generate_almanac_from_arc", &EphemerisGenApp::GenerateAlmanacFromArc,
           py::arg("t_gen_s"), py::arg("t_s"), py::arg("rv"),
           py::return_value_policy::reference_internal,
           "Fit and store an almanac message from an explicit sampled arc (`t_s` [s], "
           "`rv` [N x 6]).")
      .def("get_ephemeris_messages", &EphemerisGenApp::GetEphemerisMessages,
           py::return_value_policy::reference_internal, "All generated ephemeris messages.")
      .def("get_almanac_messages", &EphemerisGenApp::GetAlmanacMessages,
           py::return_value_policy::reference_internal, "All generated almanac messages.")
      .def("latest_ephemeris", &EphemerisGenApp::LatestEphemeris, py::arg("t_s"),
           py::return_value_policy::reference_internal,
           "Most recent ephemeris message whose validity window contains `t_s` (None if none).")
      .def("latest_almanac", &EphemerisGenApp::LatestAlmanac, py::arg("t_s"),
           py::return_value_policy::reference_internal,
           "Most recent almanac message whose validity window contains `t_s` (None if none).")
      .def("get_config", &EphemerisGenApp::GetConfig, py::return_value_policy::reference_internal,
           "The EphemerisGenConfig this app uses.");
}
