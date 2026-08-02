
#include <lupnt/lupnt.h>

#include "py_pybind11.h"
#include "py_vectorized_macros.h"

namespace py = pybind11;
using namespace lupnt;

void InitAntenna(py::module &m) {
  py::class_<Antenna>(m, "Antenna", "Antenna gain pattern loaded from a named pattern file.")
      .def(py::init<const std::string &>(), "Load the named antenna gain pattern (empty = omni).")
      .DEF_CLASS_REAL_REAL("compute_gain", Antenna, ComputeGain, "theta", "phi")
      .def(
          "compute_gain_azimuth_averaged",
          [](Antenna &ant, Real phi) { return ant.ComputeGainAzimuthAveraged(phi); },
          py::arg("phi"),
          "Azimuth-averaged gain [dB] at an off-boresight angle [rad] (linear-power "
          "average over all azimuth); NaN outside coverage, 0 for omni.")
      .def(
          "compute_gain_azimuth_averaged",
          [](Antenna &ant, const VecX &phi) { return ant.ComputeGainAzimuthAveraged(phi); },
          py::arg("phi"))
      .def(
          "get_gain_matrix", [](Antenna &ant) { return ant.GetGainMatrix().cast<double>(); },
          "Gain pattern table [dB] over (phi, theta).")
      .def(
          "get_phi_vector", [](Antenna &ant) { return ant.GetPhiVector().cast<double>(); },
          "Phi (elevation/off-boresight) sample grid [rad].")
      .def(
          "get_theta_vector", [](Antenna &ant) { return ant.GetThetaVector().cast<double>(); },
          "Theta (azimuth) sample grid [rad].")
      .def_property_readonly("name", &Antenna::GetName, "Name of the loaded antenna pattern.");
}
