
#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitTle(py::module &m) {
  py::class_<TLE>(m, "TLE", "Two-line element set describing a satellite's mean orbital elements.")
      .def_static("from_lines", &TLE::FromLines, py::arg("line1"), py::arg("line2"),
                  py::arg("line3"), "Parse a TLE from a name line plus the two element lines.")
      .def_static("from_file", &TLE::FromFile, py::arg("filename"),
                  "Parse all TLEs from a file, returning a list.")
      .def_readonly("name", &TLE::name, "Satellite name.")
      .def_readonly("prn", &TLE::prn, "PRN / catalog number.")
      .def_readonly("epoch_year", &TLE::epoch_year, "Last two digits of the epoch year [yr].")
      .def_readonly("epoch_day", &TLE::epoch_day, "Day of year and fraction of the epoch [day].")
      .def_readonly("bstar", &TLE::bstar, "B* drag term [1/R_EARTH].")
      .def_readonly("inclination", &TLE::inclination, "Inclination [deg].")
      .def_readonly("raan", &TLE::raan, "Right ascension of the ascending node [deg].")
      .def_readonly("eccentricity", &TLE::eccentricity, "Eccentricity [-].")
      .def_readonly("arg_perigee", &TLE::arg_perigee, "Argument of perigee [deg].")
      .def_readonly("mean_anomaly", &TLE::mean_anomaly, "Mean anomaly [deg].")
      .def_readonly("mean_motion", &TLE::mean_motion, "Mean motion [revs/day].")
      .def_readonly("epoch_tai", &TLE::epoch_tai, "Epoch [s TAI].")
      .def("__repr__", [](const TLE &t) { return "<TLE " + t.name + ">"; });
}
