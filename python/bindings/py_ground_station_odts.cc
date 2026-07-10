#include <lupnt/applications/ground_station/ground_station_manager_app.h>
#include <lupnt/applications/ground_station/ground_station_tracking_app.h>

#include <memory>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

// Agent-based ground-station orbit determination (Example 7). Requires the
// framework bases (Application/Agent) registered by InitSimulation first.
void InitGroundStationOdts(py::module& m) {
  // ---- Per-station sensor app (runs on each GroundStation) ----
  py::class_<GroundStationTrackingApp, Application, std::shared_ptr<GroundStationTrackingApp>>(
      m, "GroundStationTrackingApp")
      .def("station_name", &GroundStationTrackingApp::StationName)
      .def("elevation_time", &GroundStationTrackingApp::ElevationTime,
           "Sim-relative epochs [s] at which elevation was evaluated")
      .def("elevation_deg", &GroundStationTrackingApp::ElevationDeg,
           "Topocentric elevation [deg] of the target at each elevation_time")
      .def("measurement_time", &GroundStationTrackingApp::MeasurementTime,
           "Sim-relative epochs [s] of the visible (above-mask) measurements")
      .def("measurement_range", &GroundStationTrackingApp::MeasurementRange, "Noisy range [m]")
      .def("measurement_range_rate", &GroundStationTrackingApp::MeasurementRangeRate,
           "Noisy range-rate [m/s]");

  // ---- Centralized estimator app (runs on the GroundStationManager) ----
  py::class_<GroundStationManagerApp, Application, std::shared_ptr<GroundStationManagerApp>>(
      m, "GroundStationManagerApp")
      .def("has_solved", &GroundStationManagerApp::HasSolved)
      .def("converged", &GroundStationManagerApp::Converged)
      .def("num_iterations", &GroundStationManagerApp::NumIterations)
      .def("num_measurements", &GroundStationManagerApp::NumMeasurements)
      .def("station_names", &GroundStationManagerApp::StationNames)
      // Epoch-state solution + formal covariance
      .def("x0_true", &GroundStationManagerApp::X0True)
      .def("x0_initial_guess", &GroundStationManagerApp::X0InitialGuess)
      .def("x0_estimated", &GroundStationManagerApp::X0Estimated)
      .def("covariance", &GroundStationManagerApp::Covariance)
      // Full-arc time series (uniform epoch grid, world frame)
      .def("time_grid", &GroundStationManagerApp::TimeGrid, "Sim-relative epoch grid [s], size N")
      .def("truth_state", &GroundStationManagerApp::TruthState, "[N x 6] truth trajectory")
      .def("estimated_state", &GroundStationManagerApp::EstimatedState, "[N x 6] batch estimate")
      .def("estimated_covariance", &GroundStationManagerApp::EstimatedCovariance,
           "[N x 36] row-major 6x6 covariance Phi P0 Phi^T at each epoch")
      // Batch iteration history
      .def("iteration_state", &GroundStationManagerApp::IterationState)
      .def("iteration_pos_error", &GroundStationManagerApp::IterationPosError)
      .def("iteration_vel_error", &GroundStationManagerApp::IterationVelError)
      .def("iteration_correction_norm", &GroundStationManagerApp::IterationCorrectionNorm)
      .def("iteration_weighted_rms", &GroundStationManagerApp::IterationWeightedRms)
      .def("iteration_rms_range", &GroundStationManagerApp::IterationRmsRange)
      .def("iteration_rms_range_rate", &GroundStationManagerApp::IterationRmsRangeRate)
      // SRIF forward filter + smoother
      .def("srif_filtered_state", &GroundStationManagerApp::SrifFilteredState)
      .def("srif_filtered_covariance", &GroundStationManagerApp::SrifFilteredCovariance)
      .def("srif_smoothed_state", &GroundStationManagerApp::SrifSmoothedState)
      .def("srif_smoothed_covariance", &GroundStationManagerApp::SrifSmoothedCovariance);
}
