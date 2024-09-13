
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_antenna(py::module &m) {
  py::class_<Antenna>(m, "Antenna")
      .def(py::init<const std::string &>())
      .def(
          "compute_gain",
          [](Antenna &ant, double elev, double azim) -> double {
            return ant.ComputeGain(elev, azim);
          },
          py::arg("elev"), py::arg("azim"))
      .def(
          "compute_gain",
          [](Antenna &ant, Eigen::VectorXd elev, Eigen::VectorXd azim) -> Eigen::VectorXd {
            return ant.ComputeGain(elev, azim).cast<double>();
          },
          py::arg("elev"), py::arg("azim"))
      .def("get_gain_pattern", [](Antenna &ant) { return ant.GetGainPattern().cast<double>(); })
      .def("get_elevation_angles",
           [](Antenna &ant) { return ant.GetElevationAngles().cast<double>(); })
      .def("get_azimuth_angles",
           [](Antenna &ant) { return ant.GetAzimuthAngles().cast<double>(); });
}
