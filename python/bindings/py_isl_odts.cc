/**
 * @file py_isl_odts.cc
 * @brief Python bindings for the distributed ISL ODTS scenario apps.
 *
 * The scenario is built from PHYSICAL agents in `configs/isl_odts_distributed.yaml`
 * (`pnt.Simulation`): each satellite is a `Spacecraft` hosting a `SatelliteOdtsApp`
 * (onboard Schmidt-EKF) and each ground station is a `SurfaceStation` whose
 * `StationBeaconSensor` feeds a `SurfaceStationManager`'s `GroundOdtsApp` (centralized
 * ground filter). These two app types expose the per-epoch truth/estimate series read by
 * `python/examples/ex8_run_isl_odts.py`. The result accessors return `Eigen` matrices /
 * `std::vector`, converted by pybind11's built-in Eigen<->NumPy and `stl.h` list support.
 */
#include <lupnt/applications/lunar_sat_odts/ground_odts_app.h>
#include <lupnt/applications/lunar_sat_odts/satellite_odts_app.h>
#include <lupnt/lupnt.h>

#include <memory>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitIslOdts(py::module& m) {
  // ---- ISL ODTS scenario: physical `Spacecraft` onboard app + ground-segment app ----
  py::class_<SatelliteOdtsApp, Application, std::shared_ptr<SatelliteOdtsApp>>(m,
                                                                               "SatelliteOdtsApp")
      .def("sat_name", &SatelliteOdtsApp::SatName)
      .def("neighbor_names", &SatelliteOdtsApp::NeighborNames)
      .def("time_grid", &SatelliteOdtsApp::TimeGrid)
      .def("truth_state", &SatelliteOdtsApp::TruthState, "[N x 8] own truth [r,v,cb,cd]")
      .def("own_estimate", &SatelliteOdtsApp::OwnEstimate, "[N x 8] own onboard estimate")
      .def("own_cov_diag", &SatelliteOdtsApp::OwnCovDiag, "[N x 8] own covariance diagonal")
      .def("own_cov_full", &SatelliteOdtsApp::OwnCovFull, "[N x 64] own 8x8 covariance (row-major)")
      .def("num_anchors_total", &SatelliteOdtsApp::NumAnchorsTotal);

  py::class_<GroundOdtsApp, Application, std::shared_ptr<GroundOdtsApp>>(m, "GroundOdtsApp")
      .def("satellite_names", &GroundOdtsApp::SatelliteNames)
      .def("time_grid", &GroundOdtsApp::TimeGrid)
      .def("est_central", &GroundOdtsApp::EstCentral, "[N x 8*n_sat] centralized ground estimate")
      .def("cov_central_full", &GroundOdtsApp::CovCentralFull,
           "n_sat x [N x 64] per-satellite 8x8 covariance (row-major)");
}
