
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_tle(py::module &m) {
  py::class_<TLE>(m, "TLE")
      .def_static("from_lines", &TLE::FromLines, py::arg("line1"), py::arg("line2"),
                  py::arg("line3"))
      .def_static("from_file", &TLE::FromFile, py::arg("filename"))
      .def_readonly("name", &TLE::name)
      .def_readonly("prn", &TLE::prn)
      .def_readonly("epoch_year", &TLE::epoch_year)
      .def_readonly("epoch_day", &TLE::epoch_day)
      .def_readonly("bstar", &TLE::bstar)
      .def_readonly("inclination", &TLE::inclination)
      .def_readonly("raan", &TLE::raan)
      .def_readonly("eccentricity", &TLE::eccentricity)
      .def_readonly("arg_perigee", &TLE::arg_perigee)
      .def_readonly("mean_anomaly", &TLE::mean_anomaly)
      .def_readonly("mean_motion", &TLE::mean_motion)
      .def_readonly("epoch_tai", &TLE::epoch_tai)
      .def("__repr__", [](const TLE &t) { return "<TLE " + t.name + ">"; });
}
