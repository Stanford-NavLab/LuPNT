#include <lupnt/core/constants.h>
#include <lupnt/core/definitions.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitConstants(py::module& m) {
  m.def("get_lupnt_epoch", &lupnt::GetLupntEpoch);
  m.def("set_lupnt_epoch", &lupnt::SetLupntEpoch);

  py::class_<UnitSystem>(m, "UnitSystem")
      .def(py::init<double, double, double>(), py::arg("length") = METER, py::arg("time") = SECOND,
           py::arg("mass") = KILOGRAM)
      .def_readwrite("length", &UnitSystem::length)
      .def_readwrite("time", &UnitSystem::time)
      .def_readwrite("mass", &UnitSystem::mass)
      .def("from_si", &UnitSystem::FromSI, py::arg("value"), py::arg("length_power"),
           py::arg("time_power") = 0, py::arg("mass_power") = 0)
      .def("to_si", &UnitSystem::ToSI, py::arg("value"), py::arg("length_power"),
           py::arg("time_power") = 0, py::arg("mass_power") = 0)
      .def("length_from_si", &UnitSystem::Length)
      .def("area_from_si", &UnitSystem::Area)
      .def("velocity_from_si", &UnitSystem::Velocity)
      .def("acceleration_from_si", &UnitSystem::Acceleration)
      .def("gm_from_si", &UnitSystem::GravitationalParameter)
      .def("pressure_from_si", &UnitSystem::Pressure)
      .def("area_per_mass_from_si", &UnitSystem::AreaPerMass);

  py::class_<PhysicalConstants>(m, "PhysicalConstants")
      .def_readonly("GM_SUN", &PhysicalConstants::GM_SUN)
      .def_readonly("GM_EARTH", &PhysicalConstants::GM_EARTH)
      .def_readonly("GM_MOON", &PhysicalConstants::GM_MOON)
      .def_readonly("GM_MARS", &PhysicalConstants::GM_MARS)
      .def_readonly("R_SUN", &PhysicalConstants::R_SUN)
      .def_readonly("R_EARTH", &PhysicalConstants::R_EARTH)
      .def_readonly("R_MOON", &PhysicalConstants::R_MOON)
      .def_readonly("R_MARS", &PhysicalConstants::R_MARS)
      .def_readonly("WGS84_A", &PhysicalConstants::WGS84_A)
      .def_readonly("OMEGA_EARTH", &PhysicalConstants::OMEGA_EARTH)
      .def_readonly("OMEGA_MOON", &PhysicalConstants::OMEGA_MOON)
      .def_readonly("AU", &PhysicalConstants::AU)
      .def_readonly("C", &PhysicalConstants::C)
      .def_readonly("P_SUN", &PhysicalConstants::P_SUN)
      .def_readonly("coordinate_scale", &PhysicalConstants::coordinate_scale);

  m.def("get_physical_constants", py::overload_cast<const UnitSystem&>(&GetPhysicalConstants),
        py::arg("units") = SI_UNITS);
  m.def("get_physical_constants",
        py::overload_cast<const UnitSystem&, CoordinateScale>(&GetPhysicalConstants),
        py::arg("units"), py::arg("coordinate_scale"));
  m.def("get_physical_constants", py::overload_cast<CoordinateScale>(&GetPhysicalConstants),
        py::arg("coordinate_scale"));

  m.attr("METER") = py::float_(METER);
  m.attr("KILOMETER") = py::float_(KILOMETER);
  m.attr("SECOND") = py::float_(SECOND);
  m.attr("KILOGRAM") = py::float_(KILOGRAM);
  m.attr("SI_UNITS") = SI_UNITS;
  m.attr("M_S_KG_UNITS") = M_S_KG_UNITS;
  m.attr("KM_S_KG_UNITS") = KM_S_KG_UNITS;

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

  // Length
  m.attr("INCH_M") = py::float_(INCH_M);
  m.attr("FOOT_M") = py::float_(FOOT_M);
  m.attr("MILE_M") = py::float_(MILE_M);
  m.attr("KM_M") = py::float_(KM_M);
  m.attr("M_KM") = py::float_(M_KM);
  m.attr("MM_KM") = py::float_(MM_KM);
  m.attr("KM_MM") = py::float_(KM_MM);
  m.attr("MM_M") = py::float_(MM_M);
  m.attr("M_MM") = py::float_(M_MM);
  m.attr("M_CM") = py::float_(M_CM);
  m.attr("CM_M") = py::float_(CM_M);

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

  m.attr("JD_MJD_OFFSET") = py::float_(JD_MJD_OFFSET);
  m.attr("TT_TAI_OFFSET") = py::float_(TT_TAI_OFFSET);
  m.attr("A1_TAI_OFFSET") = py::float_(A1_TAI_OFFSET);

  m.attr("JD_CCSDS_TAI") = py::float_(JD_CCSDS_TAI);
  m.attr("JD_J2000_TT") = py::float_(JD_J2000_TT);
  m.attr("MJD_CCSDS_TAI") = py::float_(MJD_CCSDS_TAI);
  m.attr("MJD_J2000_TT") = py::float_(MJD_J2000_TT);

  m.attr("MJD_COORDINATE_TAI") = py::float_(MJD_COORDINATE_TAI);
  m.attr("MJD_COORDINATE_TT_TCG_TCB") = py::float_(MJD_COORDINATE_TT_TCG_TCB);

  m.attr("L_B") = py::float_(L_B);
  m.attr("L_G") = py::float_(L_G);
  m.attr("L_L") = py::float_(L_L);

  py::enum_<CoordinateScale>(m, "CoordinateScale")
      .value("TCB", CoordinateScale::TCB)
      .value("TDB", CoordinateScale::TDB)
      .value("TCG", CoordinateScale::TCG)
      .value("TT", CoordinateScale::TT)
      .value("TCL", CoordinateScale::TCL)
      .value("TL", CoordinateScale::TL)
      .export_values();

  m.def("are_coordinate_scales_convertible", &AreCoordinateScalesConvertible, py::arg("from"),
        py::arg("to"));
  m.def("coordinate_scale_factor", &CoordinateScaleFactor, py::arg("scale"));
  m.def("coordinate_scale_ratio", &CheckedCoordinateScaleRatio, py::arg("from"), py::arg("to"));
  m.def("scale_length_for_coordinate_scale", &ScaleLengthForCoordinateScaleChecked,
        py::arg("value"), py::arg("from"), py::arg("to"));
  m.def("scale_gm_for_coordinate_scale", &ScaleGravitationalParameterForCoordinateScaleChecked,
        py::arg("value"), py::arg("from"), py::arg("to"));

  // Coordinate system constants
  m.attr("GM_SUN") = py::float_(GM_SUN);
  m.attr("GM_MERCURY") = py::float_(GM_MERCURY);
  m.attr("GM_VENUS") = py::float_(GM_VENUS);
  m.attr("GM_EARTH") = py::float_(GM_EARTH);
  m.attr("GM_MOON") = py::float_(GM_MOON);
  m.attr("GM_MARS_SYSTEM") = py::float_(GM_MARS_SYSTEM);
  m.attr("GM_JUPITER_SYSTEM") = py::float_(GM_JUPITER_SYSTEM);
  m.attr("GM_SATURN_SYSTEM") = py::float_(GM_SATURN_SYSTEM);
  m.attr("GM_URANUS_SYSTEM") = py::float_(GM_URANUS_SYSTEM);
  m.attr("GM_NEPTUNE_SYSTEM") = py::float_(GM_NEPTUNE_SYSTEM);
  m.attr("GM_PLUTO_SYSTEM") = py::float_(GM_PLUTO_SYSTEM);
  m.attr("GM_CERES") = py::float_(GM_CERES);
  m.attr("GM_VESTA") = py::float_(GM_VESTA);

  // Distance
  m.attr("D_EARTH_MOON") = py::float_(D_EARTH_MOON);
  m.attr("D_EARTH_EMB") = py::float_(D_EARTH_EMB);
  m.attr("R_MOON") = py::float_(R_MOON);
  m.attr("R_EARTH") = py::float_(R_EARTH);
  m.attr("OMEGA_EARTH_MOON") = py::float_(OMEGA_EARTH_MOON);
  m.attr("D_MOON_EMB") = py::float_(D_MOON_EMB);

  m.attr("AU") = py::float_(AU);
  m.attr("C") = py::float_(C);

  py::enum_<BodyId>(m, "BodyId")
      .value("SSB", BodyId::SSB)
      .value("SOLAR_SYSTEM_BARYCENTER", BodyId::SOLAR_SYSTEM_BARYCENTER)
      .value("MERCURY_BARYCENTER", BodyId::MERCURY_BARYCENTER)
      .value("VENUS_BARYCENTER", BodyId::VENUS_BARYCENTER)
      .value("EMB", BodyId::EMB)
      .value("EARTH_MOON_BARYCENTER", BodyId::EARTH_MOON_BARYCENTER)
      .value("MARS_BARYCENTER", BodyId::MARS_BARYCENTER)
      .value("JUPITER_BARYCENTER", BodyId::JUPITER_BARYCENTER)
      .value("SATURN_BARYCENTER", BodyId::SATURN_BARYCENTER)
      .value("URANUS_BARYCENTER", BodyId::URANUS_BARYCENTER)
      .value("NEPTUNE_BARYCENTER", BodyId::NEPTUNE_BARYCENTER)
      .value("PLUTO_BARYCENTER", BodyId::PLUTO_BARYCENTER)
      .value("SUN", BodyId::SUN)
      .value("MERCURY", BodyId::MERCURY)
      .value("VENUS", BodyId::VENUS)
      .value("EARTH", BodyId::EARTH)
      .value("MOON", BodyId::MOON)
      .value("MARS", BodyId::MARS)
      .value("PHOBOS", BodyId::PHOBOS)
      .value("DEIMOS", BodyId::DEIMOS)
      .value("JUPITER", BodyId::JUPITER)
      .export_values();

  py::enum_<Time>(m, "Time")
      .value("UT1", Time::UT1)
      .value("UTC", Time::UTC)
      .value("TAI", Time::TAI)
      .value("TDB", Time::TDB)
      .value("TT", Time::TT)
      .value("TCG", Time::TCG)
      .value("TCB", Time::TCB)
      .value("GPS", Time::GPS)
      .value("JD_TT", Time::JD_TT)
      .value("JD_TDB", Time::JD_TDB)
      .value("TCL", Time::TCL)
      .value("LT", Time::LT)
      .export_values();
}
