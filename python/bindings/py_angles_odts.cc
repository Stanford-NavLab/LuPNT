/**
 * @file py_angles_odts.cc
 * @brief Python bindings for the satellite-to-satellite angles-only ODTS scenario.
 *
 * The scenario is built from PHYSICAL agents in `configs/sat_bearing_odts.yaml`
 * (`pnt.Simulation`): an OBSERVER `Spacecraft` hosting an `AnglesOdtsApp` measures the unit
 * line-of-sight (bearing) to a TARGET `Spacecraft` whose ephemeris is treated as known, and
 * runs an EKF estimating its OWN orbit 6-state `[r, v]`. The result accessors expose the
 * per-epoch truth/estimate/1-sigma series (first Monte-Carlo run) and the ensemble
 * position/velocity-error RMS, converted by pybind11's built-in Eigen<->NumPy support.
 */
#include <lupnt/applications/angles_odts/angles_odts_app.h>
#include <lupnt/lupnt.h>

#include <memory>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitAnglesOdts(py::module& m) {
  py::class_<AnglesOdtsApp, Application, std::shared_ptr<AnglesOdtsApp>>(m, "AnglesOdtsApp")
      .def("time_grid", &AnglesOdtsApp::TimeGrid)
      .def("truth_state", &AnglesOdtsApp::TruthState, "[N x 6] observer truth [r, v]")
      .def("est_state", &AnglesOdtsApp::EstState, "[N x 6] observer estimate (MC run 0)")
      .def("sigma_state", &AnglesOdtsApp::SigmaState, "[N x 6] 1-sigma (MC run 0)")
      .def("position_error_rms", &AnglesOdtsApp::PositionErrorRms, "[N] position error RMS over MC")
      .def("velocity_error_rms", &AnglesOdtsApp::VelocityErrorRms, "[N] velocity error RMS over MC")
      .def("trajectory", &AnglesOdtsApp::Trajectory,
           "[N x 19] [t, truth(6), est(6), sigma(6)] (MC run 0)")
      .def("final_position_error_m", &AnglesOdtsApp::FinalPositionErrorM)
      .def("rms_position_error_m", &AnglesOdtsApp::RmsPositionErrorM)
      .def("monte_carlo_runs", &AnglesOdtsApp::MonteCarloRuns);
}
