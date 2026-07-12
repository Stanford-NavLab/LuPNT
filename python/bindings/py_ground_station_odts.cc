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
      m, "GroundStationTrackingApp",
      "Ground-station tracking (sensor) app: generates one station's noisy two-way "
      "range/range-rate observations of a target satellite above the elevation mask and "
      "reports them to the GroundStationManagerApp")
      .def("station_name", &GroundStationTrackingApp::StationName, "Host station agent name")
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
      m, "GroundStationManagerApp",
      "Centralized ground-segment orbit-determination app: aggregates all stations' "
      "range/range-rate observations and runs a weighted-least-squares batch filter plus "
      "an SRIF forward filter and smoother over the arc")
      .def("has_solved", &GroundStationManagerApp::HasSolved,
           "True once the end-of-arc Solve has run")
      .def("converged", &GroundStationManagerApp::Converged, "True if the batch filter converged")
      .def("num_iterations", &GroundStationManagerApp::NumIterations,
           "Batch filter iteration count")
      .def("num_measurements", &GroundStationManagerApp::NumMeasurements,
           "Total aggregated observations across all stations")
      .def("station_names", &GroundStationManagerApp::StationNames,
           "Registered tracking-station names")
      // Epoch-state solution + formal covariance
      .def("x0_true", &GroundStationManagerApp::X0True, "Epoch-state truth [r,v] (world frame)")
      .def("x0_initial_guess", &GroundStationManagerApp::X0InitialGuess,
           "Perturbed a-priori epoch-state guess")
      .def("x0_estimated", &GroundStationManagerApp::X0Estimated,
           "Batch-estimated epoch state [r,v]")
      .def("covariance", &GroundStationManagerApp::Covariance,
           "6x6 formal covariance of the epoch-state estimate")
      // Full-arc time series (uniform epoch grid, world frame)
      .def("time_grid", &GroundStationManagerApp::TimeGrid, "Sim-relative epoch grid [s], size N")
      .def("truth_state", &GroundStationManagerApp::TruthState, "[N x 6] truth trajectory")
      .def("estimated_state", &GroundStationManagerApp::EstimatedState, "[N x 6] batch estimate")
      .def("estimated_covariance", &GroundStationManagerApp::EstimatedCovariance,
           "[N x 36] row-major 6x6 covariance Phi P0 Phi^T at each epoch")
      // Batch iteration history
      .def("iteration_state", &GroundStationManagerApp::IterationState,
           "[K x 6] epoch-state estimate at each batch iteration")
      .def("iteration_pos_error", &GroundStationManagerApp::IterationPosError,
           "[K] epoch-state position error [m] per iteration")
      .def("iteration_vel_error", &GroundStationManagerApp::IterationVelError,
           "[K] epoch-state velocity error [m/s] per iteration")
      .def("iteration_correction_norm", &GroundStationManagerApp::IterationCorrectionNorm,
           "[K] state-correction norm per iteration")
      .def("iteration_weighted_rms", &GroundStationManagerApp::IterationWeightedRms,
           "[K] weighted RMS of measurement residuals per iteration")
      .def("iteration_rms_range", &GroundStationManagerApp::IterationRmsRange,
           "[K] RMS range residual [m] per iteration")
      .def("iteration_rms_range_rate", &GroundStationManagerApp::IterationRmsRangeRate,
           "[K] RMS range-rate residual [m/s] per iteration")
      // SRIF forward filter + smoother
      .def("srif_filtered_state", &GroundStationManagerApp::SrifFilteredState,
           "[N x 6] SRIF forward-filtered trajectory")
      .def("srif_filtered_covariance", &GroundStationManagerApp::SrifFilteredCovariance,
           "[N x 36] row-major 6x6 SRIF forward-filtered covariance at each epoch")
      .def("srif_smoothed_state", &GroundStationManagerApp::SrifSmoothedState,
           "[N x 6] SRIF/Dyer-McReynolds smoothed trajectory")
      .def("srif_smoothed_covariance", &GroundStationManagerApp::SrifSmoothedCovariance,
           "[N x 36] row-major 6x6 smoothed covariance at each epoch");
}
