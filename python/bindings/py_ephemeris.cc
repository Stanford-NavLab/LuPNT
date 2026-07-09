/**
 * @file py_ephemeris.cc
 * @brief Python bindings for `lupnt::CartesianEphemeris`/`lupnt::Almanac`
 *        (`lupnt/applications/lunanet_ephemeris.h`, `lupnt/applications/lunanet_almanac.h`) and
 *        `lupnt::EphemerisSimulation`
 *        (`lupnt/simulations/ephemeris/ephemeris_simulation.h`).
 *
 * Config/result/options structs are plain `double`-valued (no `lupnt::Real`), so
 * they bind directly via `def_readwrite`/`def_readonly` with pybind11's built-in
 * Eigen<->NumPy conversion, following the `IslOdtsConfig`/`IslOdtsResults` pattern
 * in `py_isl_odts.cc`.
 */
#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitEphemeris(py::module& m) {
  // ---- EphemerisFitErrorStats -------------------------------------------------

  py::class_<EphemerisFitErrorStats>(m, "EphemerisFitErrorStats")
      .def(py::init<>())
      .def_readonly("rms_pos_m", &EphemerisFitErrorStats::rms_pos_m,
                    "[R, T, N, 3D] position RMS error [m]")
      .def_readonly("rms_vel_mps", &EphemerisFitErrorStats::rms_vel_mps,
                    "[R, T, N, 3D] velocity RMS error [m/s]")
      .def_readonly("p95_pos_m", &EphemerisFitErrorStats::p95_pos_m,
                    "[R, T, N, 3D] 95th-percentile |position error| [m]")
      .def_readonly("p95_vel_mps", &EphemerisFitErrorStats::p95_vel_mps,
                    "[R, T, N, 3D] 95th-percentile |velocity error| [m/s]");

  // ---- CartesianEphemeris ------------------------------------------------------

  py::class_<EphemerisFitOptions>(m, "EphemerisFitOptions")
      .def(py::init<>())
      .def_readwrite("poly_order", &EphemerisFitOptions::poly_order)
      .def_readwrite("use_keplerian_baseline", &EphemerisFitOptions::use_keplerian_baseline)
      .def_readwrite("gm", &EphemerisFitOptions::gm)
      .def_readwrite("frame", &EphemerisFitOptions::frame)
      .def_readwrite("num_fourier_terms", &EphemerisFitOptions::num_fourier_terms);

  py::class_<CartesianEphemeris>(m, "CartesianEphemeris")
      .def(py::init<>())
      .def(py::init<EphemerisFitOptions>(), py::arg("options"))
      .def("fit", &CartesianEphemeris::Fit, py::arg("t_s"), py::arg("rv"))
      .def("eval", &CartesianEphemeris::Eval, py::arg("t_s"), py::arg("params"))
      .def("eval_error", &CartesianEphemeris::EvalError, py::arg("t_s"), py::arg("rv_ref"),
           py::arg("params"))
      .def("num_params", &CartesianEphemeris::NumParams)
      .def("param_names", &CartesianEphemeris::ParamNames)
      .def("get_options", &CartesianEphemeris::GetOptions,
           py::return_value_policy::reference_internal);

  // ---- Almanac ------------------------------------------------------------------

  py::class_<AlmanacFitOptions>(m, "AlmanacFitOptions")
      .def(py::init<>())
      .def_readwrite("poly_order", &AlmanacFitOptions::poly_order)
      .def_readwrite("num_fourier_terms", &AlmanacFitOptions::num_fourier_terms)
      .def_readwrite("gm", &AlmanacFitOptions::gm)
      .def_readwrite("sidereal_period_s", &AlmanacFitOptions::sidereal_period_s)
      .def_readwrite("frame", &AlmanacFitOptions::frame);

  py::class_<Almanac>(m, "Almanac")
      .def(py::init<>())
      .def(py::init<AlmanacFitOptions>(), py::arg("options"))
      .def("fit", &Almanac::Fit, py::arg("t_s"), py::arg("rv"))
      .def("eval", &Almanac::Eval, py::arg("t_s"), py::arg("params"))
      .def("eval_error", &Almanac::EvalError, py::arg("t_s"), py::arg("rv_ref"), py::arg("params"))
      .def("num_params", &Almanac::NumParams)
      .def("param_names", &Almanac::ParamNames)
      .def("get_options", &Almanac::GetOptions, py::return_value_policy::reference_internal);

  // ---- EphemerisSimulation -------------------------------------------------------

  py::class_<EphemerisOrbitConfig>(m, "EphemerisOrbitConfig")
      .def(py::init<>())
      .def_readwrite("a_m", &EphemerisOrbitConfig::a_m)
      .def_readwrite("ecc", &EphemerisOrbitConfig::ecc)
      .def_readwrite("inc_rad", &EphemerisOrbitConfig::inc_rad)
      .def_readwrite("raan_rad", &EphemerisOrbitConfig::raan_rad)
      .def_readwrite("argp_rad", &EphemerisOrbitConfig::argp_rad)
      .def_readwrite("m0_rad", &EphemerisOrbitConfig::m0_rad)
      .def_readwrite("coe_frame", &EphemerisOrbitConfig::coe_frame);

  py::class_<EphemerisSimulationConfig>(m, "EphemerisSimulationConfig")
      .def(py::init<>())
      .def_readwrite("start_epoch_utc", &EphemerisSimulationConfig::start_epoch_utc)
      .def_readwrite("orbit", &EphemerisSimulationConfig::orbit)
      .def_readwrite("propagate_frame", &EphemerisSimulationConfig::propagate_frame)
      .def_readwrite("duration_days", &EphemerisSimulationConfig::duration_days)
      .def_readwrite("sample_dt_s", &EphemerisSimulationConfig::sample_dt_s)
      .def_readwrite("moon_gravity_degree", &EphemerisSimulationConfig::moon_gravity_degree)
      .def_readwrite("moon_gravity_order", &EphemerisSimulationConfig::moon_gravity_order)
      .def_readwrite("include_earth", &EphemerisSimulationConfig::include_earth)
      .def_readwrite("include_sun", &EphemerisSimulationConfig::include_sun)
      .def_readwrite("use_relativity", &EphemerisSimulationConfig::use_relativity)
      .def_readwrite("integration_step_s", &EphemerisSimulationConfig::integration_step_s)
      .def_readwrite("fit_window_minutes", &EphemerisSimulationConfig::fit_window_minutes)
      .def_readwrite("num_windows", &EphemerisSimulationConfig::num_windows)
      .def_readwrite("cartesian_poly_order", &EphemerisSimulationConfig::cartesian_poly_order)
      .def_readwrite("cartesian_use_keplerian_baseline",
                     &EphemerisSimulationConfig::cartesian_use_keplerian_baseline)
      .def_readwrite("cartesian_num_fourier_terms",
                     &EphemerisSimulationConfig::cartesian_num_fourier_terms)
      .def_readwrite("almanac_poly_order", &EphemerisSimulationConfig::almanac_poly_order)
      .def_readwrite("almanac_num_fourier_terms",
                     &EphemerisSimulationConfig::almanac_num_fourier_terms)
      .def_readwrite("output_frame", &EphemerisSimulationConfig::output_frame)
      .def_readwrite("datasize_precision_m", &EphemerisSimulationConfig::datasize_precision_m);

  py::class_<EphemerisWindowResult>(m, "EphemerisWindowResult")
      .def(py::init<>())
      .def_readonly("fit_window_min", &EphemerisWindowResult::fit_window_min)
      .def_readonly("num_params", &EphemerisWindowResult::num_params)
      .def_readonly("total_bits", &EphemerisWindowResult::total_bits)
      .def_readonly("pos_rms_m", &EphemerisWindowResult::pos_rms_m)
      .def_readonly("vel_rms_mps", &EphemerisWindowResult::vel_rms_mps)
      .def_readonly("pos_p95_m", &EphemerisWindowResult::pos_p95_m)
      .def_readonly("vel_p95_mps", &EphemerisWindowResult::vel_p95_mps);

  py::class_<EphemerisSimulation>(m, "EphemerisSimulation")
      .def(py::init<EphemerisSimulationConfig>(), py::arg("config"))
      .def("setup", &EphemerisSimulation::Setup)
      .def("run", &EphemerisSimulation::Run)
      .def("get_config", &EphemerisSimulation::GetConfig,
           py::return_value_policy::reference_internal)
      .def("get_truth_times", &EphemerisSimulation::GetTruthTimes,
           py::return_value_policy::reference_internal)
      .def("get_truth_states", &EphemerisSimulation::GetTruthStates,
           py::return_value_policy::reference_internal)
      .def("get_cartesian_results", &EphemerisSimulation::GetCartesianResults,
           py::return_value_policy::reference_internal)
      .def("get_almanac_results", &EphemerisSimulation::GetAlmanacResults,
           py::return_value_policy::reference_internal);

  // ---- EphemerisGenApp (LunaNet nav-message generation sub-app) ------------------

  py::class_<BroadcastMessage>(m, "BroadcastMessage")
      .def(py::init<>())
      .def_readonly("t_generated_s", &BroadcastMessage::t_generated_s)
      .def_readonly("t_start_s", &BroadcastMessage::t_start_s)
      .def_readonly("t_end_s", &BroadcastMessage::t_end_s)
      .def_readonly("frame", &BroadcastMessage::frame)
      .def_readonly("params", &BroadcastMessage::params);

  py::class_<EphemerisGenConfig>(m, "EphemerisGenConfig")
      .def(py::init<>())
      .def_readwrite("generate_ephemeris", &EphemerisGenConfig::generate_ephemeris)
      .def_readwrite("generate_almanac", &EphemerisGenConfig::generate_almanac)
      .def_readwrite("ephemeris_options", &EphemerisGenConfig::ephemeris_options)
      .def_readwrite("ephemeris_window_s", &EphemerisGenConfig::ephemeris_window_s)
      .def_readwrite("ephemeris_refresh_s", &EphemerisGenConfig::ephemeris_refresh_s)
      .def_readwrite("almanac_options", &EphemerisGenConfig::almanac_options)
      .def_readwrite("almanac_window_s", &EphemerisGenConfig::almanac_window_s)
      .def_readwrite("almanac_refresh_s", &EphemerisGenConfig::almanac_refresh_s)
      .def_readwrite("ephemeris_fit_samples", &EphemerisGenConfig::ephemeris_fit_samples)
      .def_readwrite("almanac_fit_samples", &EphemerisGenConfig::almanac_fit_samples)
      .def_readwrite("output_frame", &EphemerisGenConfig::output_frame);

  py::class_<EphemerisGenApp>(m, "EphemerisGenApp")
      .def(py::init<>())
      .def(py::init<EphemerisGenConfig>(), py::arg("config"))
      .def("generate_ephemeris_from_arc", &EphemerisGenApp::GenerateEphemerisFromArc,
           py::arg("t_gen_s"), py::arg("t_s"), py::arg("rv"),
           py::return_value_policy::reference_internal)
      .def("generate_almanac_from_arc", &EphemerisGenApp::GenerateAlmanacFromArc,
           py::arg("t_gen_s"), py::arg("t_s"), py::arg("rv"),
           py::return_value_policy::reference_internal)
      .def("get_ephemeris_messages", &EphemerisGenApp::GetEphemerisMessages,
           py::return_value_policy::reference_internal)
      .def("get_almanac_messages", &EphemerisGenApp::GetAlmanacMessages,
           py::return_value_policy::reference_internal)
      .def("latest_ephemeris", &EphemerisGenApp::LatestEphemeris, py::arg("t_s"),
           py::return_value_policy::reference_internal)
      .def("latest_almanac", &EphemerisGenApp::LatestAlmanac, py::arg("t_s"),
           py::return_value_policy::reference_internal)
      .def("get_config", &EphemerisGenApp::GetConfig, py::return_value_policy::reference_internal);
}
