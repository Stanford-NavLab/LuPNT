#include <lupnt/core/constants.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_constants(py::module &m) {
  m.attr("PI") = py::float_(PI);
  m.attr("TWO_PI") = py::float_(TWO_PI);
  m.attr("PI_OVER_TWO") = py::float_(PI_OVER_TWO);
  m.attr("E") = py::float_(E);
  m.attr("RAD_PER_DEG") = py::float_(RAD);
  m.attr("DEG_PER_RAD") = py::float_(DEG);

  // Time System Constants
  m.attr("SECS_PER_DAY") = py::float_(SECS_PER_DAY);
  m.attr("SECS_PER_HOUR") = py::float_(SECS_HOUR);
  m.attr("SECS_PER_MINUTE") = py::float_(SECS_MINUTE);
  m.attr("DAYS_PER_YEAR") = py::float_(DAYS_YEAR);
  m.attr("DAYS_PER_JULIAN_CENTURY") = py::float_(DAYS_JULIAN_CENTURY);
  m.attr("DAYS_PER_SEC") = py::float_(DAYS_PER_SEC);
  m.attr("TIME_OF_J2000") = py::float_(TIME_OF_J2000);
  m.attr("JD_OF_J2000") = py::float_(JD_OF_J2000);
  m.attr("MJD_OF_J2000") = py::float_(MJD_J2000);
  m.attr("A1MJD_OF_J2000") = py::float_(A1MJD_OF_J2000);
  m.attr("JD_MJD_OFFSET") = py::float_(JD_MJD_OFFSET);
  m.attr("TT_TAI_OFFSET") = py::float_(TT_TAI_OFFSET);
  m.attr("A1_TAI_OFFSET") = py::float_(A1_TAI_OFFSET);
  m.attr("JD_NOV_17_1858") = py::float_(JD_NOV_17_1858);
  m.attr("JD_JAN_5_1941") = py::float_(JD_JAN_5_1941);

  // Coordinate System Constants
  m.attr("d_E_M") = py::float_(d_E_M);
  m.attr("MU_EARTH") = py::float_(MU_EARTH);
  m.attr("MU_MOON") = py::float_(GM_MOON);
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

  py::enum_<NaifId>(m, "NaifId")
      .value("SOLAR_SYSTEM_BARYCENTER", NaifId::SOLAR_SYSTEM_BARYCENTER)
      .value("SSB", NaifId::SSB)
      .value("MERCURY_BARYCENTER", NaifId::MERCURY_BARYCENTER)
      .value("VENUS_BARYCENTER", NaifId::VENUS_BARYCENTER)
      .value("EARTH_BARYCENTER", NaifId::EARTH_BARYCENTER)
      .value("EARTH_MOON_BARYCENTER", NaifId::EARTH_MOON_BARYCENTER)
      .value("MARS_BARYCENTER", NaifId::MARS_BARYCENTER)
      .value("JUPITER_BARYCENTER", NaifId::JUPITER_BARYCENTER)
      .value("SATURN_BARYCENTER", NaifId::SATURN_BARYCENTER)
      .value("URANUS_BARYCENTER", NaifId::URANUS_BARYCENTER)
      .value("NEPTUNE_BARYCENTER", NaifId::NEPTUNE_BARYCENTER)
      .value("PLUTO_BARYCENTER", NaifId::PLUTO_BARYCENTER)
      .value("SUN", NaifId::SUN)
      .value("MERCURY", NaifId::MERCURY)
      .value("VENUS", NaifId::VENUS)
      .value("EARTH", NaifId::EARTH)
      .value("MOON", NaifId::MOON)
      .value("MARS", NaifId::MARS)
      .value("PHOBOS", NaifId::PHOBOS)
      .value("DEIMOS", NaifId::DEIMOS)
      .value("JUPITER", NaifId::JUPITER)
      .export_values();
}