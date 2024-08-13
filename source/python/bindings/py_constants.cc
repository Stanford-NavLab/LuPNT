#include <lupnt/core/constants.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_constants(py::module& m) {
  // Math constants
  m.attr("PI") = py::float_(PI);
  m.attr("TWO_PI") = py::float_(TWO_PI);
  m.attr("PI_OVER_TWO") = py::float_(PI_OVER_TWO);
  m.attr("E") = py::float_(E);
  m.attr("EPS") = py::float_(EPS);

  // Angle conversion
  m.attr("RAD") = py::float_(RAD);
  m.attr("DEG") = py::float_(DEG);
  m.attr("ARCSEC_DEG") = py::float_(ARCSEC_DEG);
  m.attr("DEG_ARCSEC") = py::float_(DEG_ARCSEC);
  m.attr("RAD_ARCSEC") = py::float_(RAD_ARCSEC);
  m.attr("ARCSEC_RAD") = py::float_(ARCSEC_RAD);

  // Mass
  m.attr("LBM_TO_KG") = py::float_(LBM_TO_KG);
  m.attr("SLUG_TO_KG") = py::float_(SLUG_TO_KG);

  // Length
  m.attr("INCH_M") = py::float_(INCH_M);
  m.attr("FOOT_M") = py::float_(FOOT_M);
  m.attr("MILE_M") = py::float_(MILE_M);
  m.attr("KM_M") = py::float_(KM_M);
  m.attr("M_KM") = py::float_(M_KM);

  // Time system constants
  m.attr("SECS_DAY") = py::float_(SECS_DAY);
  m.attr("SECS_HOUR") = py::float_(SECS_HOUR);
  m.attr("SECS_MINUTE") = py::float_(SECS_MINUTE);
  m.attr("MINS_HOUR") = py::float_(MINS_HOUR);
  m.attr("MINS_DAY") = py::float_(MINS_DAY);
  m.attr("HOURS_DAY") = py::float_(HOURS_DAY);
  m.attr("DAYS_WEEK") = py::float_(DAYS_WEEK);
  m.attr("DAYS_YEAR") = py::float_(DAYS_YEAR);
  m.attr("DAYS_CENTURY") = py::float_(DAYS_CENTURY);
  m.attr("DAYS_SEC") = py::float_(DAYS_SEC);
  m.attr("TIME_OF_J2000") = py::float_(TIME_OF_J2000);
  m.attr("JD_J2000") = py::float_(JD_J2000);
  m.attr("MJD_J2000") = py::float_(MJD_J2000);
  m.attr("JD_T0") = py::float_(JD_T0);
  m.attr("JD_MJD_OFFSET") = py::float_(JD_MJD_OFFSET);
  m.attr("TT_TAI_OFFSET") = py::float_(TT_TAI_OFFSET);
  m.attr("A1_TAI_OFFSET") = py::float_(A1_TAI_OFFSET);
  m.attr("JD_JAN_5_1941") = py::float_(JD_JAN_5_1941);
  m.attr("JD_NOV_17_1858") = py::float_(JD_NOV_17_1858);
  m.attr("L_B") = py::float_(L_B);
  m.attr("L_G") = py::float_(L_G);
  m.attr("NUM_SECS") = py::float_(NUM_SECS);
  m.attr("JULIAN_DATE_OF_010541") = py::int_(JULIAN_DATE_OF_010541);

  // Coordinate system constants
  m.attr("GM_SUN") = py::float_(GM_SUN);
  m.attr("GM_MERCURY") = py::float_(GM_MERCURY);
  m.attr("GM_VENUS") = py::float_(GM_VENUS);
  m.attr("GM_EARTH") = py::float_(GM_EARTH);
  m.attr("GM_MOON") = py::float_(GM_MOON);

  m.attr("R_MOON") = py::float_(R_MOON);
  m.attr("R_EARTH") = py::float_(R_EARTH);

  py::enum_<NaifId>(m, "NaifId")
      .value("SOLAR_SYSTEM_BARYCENTER", NaifId::SOLAR_SYSTEM_BARYCENTER)
      .value("SSB", NaifId::SSB)
      .value("MERCURY_BARYCENTER", NaifId::MERCURY_BARYCENTER)
      .value("VENUS_BARYCENTER", NaifId::VENUS_BARYCENTER)
      .value("EMB", NaifId::EMB)
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

  m.attr("UT1") = py::str(TimeSys::UT1);
  m.attr("UTC") = py::str(TimeSys::UTC);
  m.attr("TAI") = py::str(TimeSys::TAI);
  m.attr("TDB") = py::str(TimeSys::TDB);
  m.attr("TT") = py::str(TimeSys::TT);
  m.attr("TCG") = py::str(TimeSys::TCG);
  m.attr("TCB") = py::str(TimeSys::TCB);
  m.attr("GPS") = py::str(TimeSys::GPS);
  m.attr("JD_TT") = py::str(TimeSys::JD_TT);
  m.attr("JD_TDB") = py::str(TimeSys::JD_TDB);
}
