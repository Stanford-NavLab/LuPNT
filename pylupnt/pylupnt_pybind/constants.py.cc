#include <lupnt/core/constants.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_constants(py::module &m) {
  m.attr("PI") = py::float_(PI);
  m.attr("TWO_PI") = py::float_(TWO_PI);
  m.attr("PI_OVER_TWO") = py::float_(PI_OVER_TWO);
  m.attr("E") = py::float_(E);
  m.attr("RAD_PER_DEG") = py::float_(RAD_PER_DEG);
  m.attr("DEG_PER_RAD") = py::float_(DEG_PER_RAD);

  // Time System Constants
  m.attr("SECS_PER_DAY") = py::float_(SECS_PER_DAY);
  m.attr("SECS_PER_HOUR") = py::float_(SECS_PER_HOUR);
  m.attr("SECS_PER_MIN") = py::float_(SECS_PER_MINUTE);
  m.attr("DAYS_PER_YEAR") = py::float_(DAYS_PER_YEAR);
  m.attr("DAYS_PER_JULIAN_CENTURY") = py::float_(DAYS_PER_JULIAN_CENTURY);
  m.attr("DAYS_PER_SEC") = py::float_(DAYS_PER_SEC);
  m.attr("TIME_OF_J2000") = py::float_(TIME_OF_J2000);
  m.attr("JD_OF_J2000") = py::float_(JD_OF_J2000);
  m.attr("MJD_OF_J2000") = py::float_(MJD_OF_J2000);
  m.attr("A1MJD_OF_J2000") = py::float_(A1MJD_OF_J2000);
  m.attr("JD_MJD_OFFSET") = py::float_(JD_MJD_OFFSET);
  m.attr("TT_TAI_OFFSET") = py::float_(TT_TAI_OFFSET);
  m.attr("A1_TAI_OFFSET") = py::float_(A1_TAI_OFFSET);
  m.attr("JD_NOV_17_1858") = py::float_(JD_NOV_17_1858);
  m.attr("JD_JAN_5_1941") = py::float_(JD_JAN_5_1941);

  // Coordinate System Constants
  m.attr("d_E_M") = py::float_(d_E_M);
  m.attr("MU_EARTH") = py::float_(MU_EARTH);
  m.attr("MU_MOON") = py::float_(MU_MOON);
  m.attr("d_E_EMB") = py::float_(d_E_EMB);
  m.attr("R_EARTH") = py::float_(R_EARTH);
  m.attr("R_MOON") = py::float_(R_MOON);
  m.attr("OMEGA_E_M") = py::float_(OMEGA_E_M);
  m.attr("d_M_EMB") = py::float_(d_M_EMB);

  m.attr("J2_EARTH") = py::float_(J2_EARTH);
  m.attr("J2_MOON") = py::float_(J2_MOON);
  m.attr("C22_MOON") = py::float_(C22_MOON);

  // Solar radiation pressure constants
  m.attr("AU") = py::float_(AU);
  m.attr("S_AU") = py::float_(S_AU);
  m.attr("C") = py::float_(C);
  m.attr("P_SUN") = py::float_(P_SUN);
}