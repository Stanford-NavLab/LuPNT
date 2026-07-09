#include <lupnt/simulations/ground_station_odts/ground_station_odts_simulation.h>

#include "py_pybind11.h"

void InitGroundStationOdts(py::module& m) {
  py::class_<GroundStationOdtsStationConfig>(m, "GroundStationOdtsStationConfig")
      .def(py::init<>())
      .def(py::init(
               [](std::string name, double latitude_deg, double longitude_deg, double altitude_m) {
                 GroundStationOdtsStationConfig cfg;
                 cfg.name = std::move(name);
                 cfg.latitude_deg = latitude_deg;
                 cfg.longitude_deg = longitude_deg;
                 cfg.altitude_m = altitude_m;
                 return cfg;
               }),
           py::arg("name"), py::arg("latitude_deg"), py::arg("longitude_deg"),
           py::arg("altitude_m"))
      .def_readwrite("name", &GroundStationOdtsStationConfig::name)
      .def_readwrite("latitude_deg", &GroundStationOdtsStationConfig::latitude_deg)
      .def_readwrite("longitude_deg", &GroundStationOdtsStationConfig::longitude_deg)
      .def_readwrite("altitude_m", &GroundStationOdtsStationConfig::altitude_m);

  m.def("default_dsn_stations", &DefaultDsnStations,
        "Default 3-station DSN 70 m subnet: Goldstone (DSS-14), Canberra (DSS-43), "
        "Madrid (DSS-63)");

  py::class_<GroundStationOdtsConfig>(m, "GroundStationOdtsConfig")
      .def(py::init<>())
      .def_readwrite("seed", &GroundStationOdtsConfig::seed)
      .def_readwrite("start_epoch_utc", &GroundStationOdtsConfig::start_epoch_utc)
      .def_readwrite("orbit_a_m", &GroundStationOdtsConfig::orbit_a_m)
      .def_readwrite("orbit_ecc", &GroundStationOdtsConfig::orbit_ecc)
      .def_readwrite("orbit_inc_rad", &GroundStationOdtsConfig::orbit_inc_rad)
      .def_readwrite("orbit_raan_rad", &GroundStationOdtsConfig::orbit_raan_rad)
      .def_readwrite("orbit_argp_rad", &GroundStationOdtsConfig::orbit_argp_rad)
      .def_readwrite("orbit_mean_anomaly_rad", &GroundStationOdtsConfig::orbit_mean_anomaly_rad)
      .def_readwrite("moon_gravity_degree", &GroundStationOdtsConfig::moon_gravity_degree)
      .def_readwrite("moon_gravity_order", &GroundStationOdtsConfig::moon_gravity_order)
      .def_readwrite("include_earth", &GroundStationOdtsConfig::include_earth)
      .def_readwrite("include_sun", &GroundStationOdtsConfig::include_sun)
      .def_readwrite("use_relativity", &GroundStationOdtsConfig::use_relativity)
      .def_readwrite("use_srp", &GroundStationOdtsConfig::use_srp)
      .def_readwrite("srp_cr", &GroundStationOdtsConfig::srp_cr)
      .def_readwrite("srp_area_m2", &GroundStationOdtsConfig::srp_area_m2)
      .def_readwrite("srp_mass_kg", &GroundStationOdtsConfig::srp_mass_kg)
      .def_readwrite("integration_step_s", &GroundStationOdtsConfig::integration_step_s)
      .def_readwrite("ground_stations", &GroundStationOdtsConfig::ground_stations)
      .def_readwrite("elevation_mask_deg", &GroundStationOdtsConfig::elevation_mask_deg)
      .def_readwrite("duration_s", &GroundStationOdtsConfig::duration_s)
      .def_readwrite("obs_interval_s", &GroundStationOdtsConfig::obs_interval_s)
      .def_readwrite("use_range", &GroundStationOdtsConfig::use_range)
      .def_readwrite("use_range_rate", &GroundStationOdtsConfig::use_range_rate)
      .def_readwrite("range_sigma_m", &GroundStationOdtsConfig::range_sigma_m)
      .def_readwrite("range_rate_sigma_mps", &GroundStationOdtsConfig::range_rate_sigma_mps)
      .def_readwrite("initial_position_sigma_m", &GroundStationOdtsConfig::initial_position_sigma_m)
      .def_readwrite("initial_velocity_sigma_mps",
                     &GroundStationOdtsConfig::initial_velocity_sigma_mps)
      .def_readwrite("batch_max_iterations", &GroundStationOdtsConfig::batch_max_iterations)
      .def_readwrite("batch_convergence_tol", &GroundStationOdtsConfig::batch_convergence_tol)
      .def_readwrite("batch_use_weights", &GroundStationOdtsConfig::batch_use_weights)
      .def_readwrite("batch_use_initialization", &GroundStationOdtsConfig::batch_use_initialization)
      .def_readwrite("batch_use_analytic_jacobian",
                     &GroundStationOdtsConfig::batch_use_analytic_jacobian)
      .def_readwrite("run_srif", &GroundStationOdtsConfig::run_srif)
      .def_readwrite("srif_use_process_noise", &GroundStationOdtsConfig::srif_use_process_noise)
      .def_readwrite("srif_accel_psd", &GroundStationOdtsConfig::srif_accel_psd);

  py::class_<GroundStationOdtsResults>(m, "GroundStationOdtsResults")
      .def(py::init<>())
      .def_readonly("station_names", &GroundStationOdtsResults::station_names)
      .def_readonly("t_tdb", &GroundStationOdtsResults::t_tdb)
      .def_readonly("truth_state", &GroundStationOdtsResults::truth_state)
      .def_readonly("elevation_deg", &GroundStationOdtsResults::elevation_deg)
      .def_readonly("obs_epoch_index", &GroundStationOdtsResults::obs_epoch_index)
      .def_readonly("obs_station_index", &GroundStationOdtsResults::obs_station_index)
      .def_readonly("obs_range_true_m", &GroundStationOdtsResults::obs_range_true_m)
      .def_readonly("obs_range_rate_true_mps", &GroundStationOdtsResults::obs_range_rate_true_mps)
      .def_readonly("obs_range_m", &GroundStationOdtsResults::obs_range_m)
      .def_readonly("obs_range_rate_mps", &GroundStationOdtsResults::obs_range_rate_mps)
      .def_readonly("iteration_state_estimate", &GroundStationOdtsResults::iteration_state_estimate)
      .def_readonly("iteration_correction_norm",
                    &GroundStationOdtsResults::iteration_correction_norm)
      .def_readonly("iteration_weighted_rms", &GroundStationOdtsResults::iteration_weighted_rms)
      .def_readonly("iteration_pos_error_m", &GroundStationOdtsResults::iteration_pos_error_m)
      .def_readonly("iteration_vel_error_mps", &GroundStationOdtsResults::iteration_vel_error_mps)
      .def_readonly("iteration_rms_range_m", &GroundStationOdtsResults::iteration_rms_range_m)
      .def_readonly("iteration_rms_range_rate_mps",
                    &GroundStationOdtsResults::iteration_rms_range_rate_mps)
      .def_readonly("x0_true", &GroundStationOdtsResults::x0_true)
      .def_readonly("x0_initial_guess", &GroundStationOdtsResults::x0_initial_guess)
      .def_readonly("x0_estimated", &GroundStationOdtsResults::x0_estimated)
      .def_readonly("covariance", &GroundStationOdtsResults::covariance)
      .def_readonly("converged", &GroundStationOdtsResults::converged)
      .def_readonly("num_iterations", &GroundStationOdtsResults::num_iterations)
      .def_readonly("estimated_state", &GroundStationOdtsResults::estimated_state)
      .def_readonly("estimated_covariance", &GroundStationOdtsResults::estimated_covariance)
      .def_readonly("srif_filtered_state", &GroundStationOdtsResults::srif_filtered_state)
      .def_readonly("srif_filtered_covariance", &GroundStationOdtsResults::srif_filtered_covariance)
      .def_readonly("srif_smoothed_state", &GroundStationOdtsResults::srif_smoothed_state)
      .def_readonly("srif_smoothed_covariance",
                    &GroundStationOdtsResults::srif_smoothed_covariance);

  py::class_<GroundStationOdtsSimulation>(m, "GroundStationOdtsSimulation")
      .def(py::init<GroundStationOdtsConfig>(), py::arg("config"))
      .def("setup", &GroundStationOdtsSimulation::Setup)
      .def("precompute", &GroundStationOdtsSimulation::Precompute)
      .def("run", &GroundStationOdtsSimulation::Run)
      .def("get_config", &GroundStationOdtsSimulation::GetConfig,
           py::return_value_policy::reference_internal)
      .def("get_results", &GroundStationOdtsSimulation::GetResults,
           py::return_value_policy::reference_internal);
}
