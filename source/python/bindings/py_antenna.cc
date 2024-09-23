
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "py_vectorized_macros.cc"

namespace py = pybind11;
using namespace lupnt;

void init_antenna(py::module &m) {
  py::class_<Antenna>(m, "Antenna")
      .def(py::init<const std::string &>())
      .DEF_CLASS_REAL_REAL("compute_gain", Antenna, ComputeGain, "theta", "phi")
      .def("get_gain_matrix", [](Antenna &ant) { return ant.GetGainMatrix().cast<double>(); })
      .def("get_phi_vector", [](Antenna &ant) { return ant.GetPhiVector().cast<double>(); })
      .def("get_theta_vector", [](Antenna &ant) { return ant.GetThetaVector().cast<double>(); })
      .def_property_readonly("name", &Antenna::GetName);
}
