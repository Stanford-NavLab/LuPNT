#include <lupnt/core/constants.h>
#include <lupnt/core/definitions.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitConstants(py::module& m) {
  m.def("get_lupnt_epoch", &lupnt::GetLupntEpoch,
        "Global LuPNT reference epoch [s past J2000, TDB] that propagation time is added to.");
  m.def("set_lupnt_epoch", &lupnt::SetLupntEpoch,
        "Set the global LuPNT reference epoch [s past J2000, TDB].");

  py::class_<UnitSystem>(m, "UnitSystem",
                         "Length/time/mass scale factors defining a working unit system.")
      .def(py::init<double, double, double>(), py::arg("length") = METER, py::arg("time") = SECOND,
           py::arg("mass") = KILOGRAM)
      .def_readwrite("length", &UnitSystem::length, "Length unit as a multiple of the SI meter.")
      .def_readwrite("time", &UnitSystem::time, "Time unit as a multiple of the SI second.")
      .def_readwrite("mass", &UnitSystem::mass, "Mass unit as a multiple of the SI kilogram.")
      .def("from_si", &UnitSystem::FromSI, py::arg("value"), py::arg("length_power"),
           py::arg("time_power") = 0, py::arg("mass_power") = 0,
           "Convert an SI value into this unit system given length/time/mass dimensional powers.")
      .def("to_si", &UnitSystem::ToSI, py::arg("value"), py::arg("length_power"),
           py::arg("time_power") = 0, py::arg("mass_power") = 0,
           "Convert a value in this unit system back to SI (inverse of from_si).")
      .def("length_from_si", &UnitSystem::Length, "Convert a length [m] into this unit system.")
      .def("area_from_si", &UnitSystem::Area, "Convert an area [m^2] into this unit system.")
      .def("velocity_from_si", &UnitSystem::Velocity,
           "Convert a velocity [m/s] into this unit system.")
      .def("acceleration_from_si", &UnitSystem::Acceleration,
           "Convert an acceleration [m/s^2] into this unit system.")
      .def("gm_from_si", &UnitSystem::GravitationalParameter,
           "Convert a gravitational parameter [m^3/s^2] into this unit system.")
      .def("pressure_from_si", &UnitSystem::Pressure,
           "Convert a pressure [Pa] into this unit system.")
      .def("area_per_mass_from_si", &UnitSystem::AreaPerMass,
           "Convert an area-per-mass [m^2/kg] into this unit system.");

  py::class_<PhysicalConstants>(m, "PhysicalConstants",
                                "Physical constants expressed in a chosen unit system.")
      .def_readonly("GM_SUN", &PhysicalConstants::GM_SUN, "Sun gravitational parameter.")
      .def_readonly("GM_EARTH", &PhysicalConstants::GM_EARTH, "Earth gravitational parameter.")
      .def_readonly("GM_MOON", &PhysicalConstants::GM_MOON, "Moon gravitational parameter.")
      .def_readonly("GM_MARS", &PhysicalConstants::GM_MARS, "Mars gravitational parameter.")
      .def_readonly("R_SUN", &PhysicalConstants::R_SUN, "Sun mean radius.")
      .def_readonly("R_EARTH", &PhysicalConstants::R_EARTH, "Earth mean radius.")
      .def_readonly("R_MOON", &PhysicalConstants::R_MOON, "Moon mean radius.")
      .def_readonly("R_MARS", &PhysicalConstants::R_MARS, "Mars mean radius.")
      .def_readonly("WGS84_A", &PhysicalConstants::WGS84_A, "WGS84 Earth equatorial radius.")
      .def_readonly("OMEGA_EARTH", &PhysicalConstants::OMEGA_EARTH,
                    "Earth rotation rate [rad/time].")
      .def_readonly("OMEGA_MOON", &PhysicalConstants::OMEGA_MOON, "Moon rotation rate [rad/time].")
      .def_readonly("AU", &PhysicalConstants::AU, "Astronomical unit.")
      .def_readonly("C", &PhysicalConstants::C, "Speed of light.")
      .def_readonly("P_SUN", &PhysicalConstants::P_SUN, "Solar radiation pressure at 1 AU.")
      .def_readonly("coordinate_scale", &PhysicalConstants::coordinate_scale,
                    "Coordinate time scale these constants are expressed in.");

  m.def("get_physical_constants", py::overload_cast<const UnitSystem&>(&GetPhysicalConstants),
        py::arg("units") = SI_UNITS, "Physical constants rescaled into the given unit system.");
  m.def("get_physical_constants",
        py::overload_cast<const UnitSystem&, CoordinateScale>(&GetPhysicalConstants),
        py::arg("units"), py::arg("coordinate_scale"),
        "Physical constants in the given unit system and coordinate time scale.");
  m.def("get_physical_constants", py::overload_cast<CoordinateScale>(&GetPhysicalConstants),
        py::arg("coordinate_scale"),
        "SI physical constants rescaled from TDB to the given coordinate time scale.");

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

  py::enum_<CoordinateScale>(m, "CoordinateScale",
                             "Relativistic coordinate time scale of a quantity.")
      .value("TCB", CoordinateScale::TCB, "Barycentric coordinate time (TCB).")
      .value("TDB", CoordinateScale::TDB, "Barycentric dynamical time (TDB).")
      .value("TCG", CoordinateScale::TCG, "Geocentric coordinate time (TCG).")
      .value("TT", CoordinateScale::TT, "Terrestrial time (TT).")
      .value("TCL", CoordinateScale::TCL, "Lunicentric coordinate time (TCL).")
      .value("TL", CoordinateScale::TL, "Lunar time (TL).")
      .export_values();

  m.def("are_coordinate_scales_convertible", &AreCoordinateScalesConvertible, py::arg("from"),
        py::arg("to"), "True if the two coordinate scales are related by a constant factor.");
  m.def("coordinate_scale_factor", &CoordinateScaleFactor, py::arg("scale"),
        "Rate factor of a coordinate scale relative to its proper coordinate time.");
  m.def(
      "coordinate_scale_ratio", &CheckedCoordinateScaleRatio, py::arg("from"), py::arg("to"),
      "Ratio to convert a quantity from one coordinate scale to another; throws if incompatible.");
  m.def("scale_length_for_coordinate_scale", &ScaleLengthForCoordinateScaleChecked,
        py::arg("value"), py::arg("from"), py::arg("to"),
        "Rescale a length from one coordinate scale to another; throws if incompatible.");
  m.def("scale_gm_for_coordinate_scale", &ScaleGravitationalParameterForCoordinateScaleChecked,
        py::arg("value"), py::arg("from"), py::arg("to"),
        "Rescale a gravitational parameter from one coordinate scale to another.");

  // Coordinate system constants
  m.attr("GM_SUN") = py::float_(GM_SUN);
  m.attr("GM_ASTEROID_BELT") = GM_ASTEROID_BELT;
  m.attr("GM_KUIPER_BELT") = GM_KUIPER_BELT;
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
  m.attr("SOLAR_FLUX_AU") = py::float_(SOLAR_FLUX_AU);
  m.attr("P_SUN") = py::float_(P_SUN);

  py::enum_<BodyId>(m, "BodyId", "NAIF-style identifier for a solar-system body or barycenter.")
      .value("SSB", BodyId::SSB, "Solar system barycenter.")
      .value("SOLAR_SYSTEM_BARYCENTER", BodyId::SOLAR_SYSTEM_BARYCENTER, "Solar system barycenter.")
      .value("MERCURY_BARYCENTER", BodyId::MERCURY_BARYCENTER, "Mercury barycenter.")
      .value("VENUS_BARYCENTER", BodyId::VENUS_BARYCENTER, "Venus barycenter.")
      .value("EMB", BodyId::EMB, "Earth-Moon barycenter.")
      .value("EARTH_MOON_BARYCENTER", BodyId::EARTH_MOON_BARYCENTER, "Earth-Moon barycenter.")
      .value("MARS_BARYCENTER", BodyId::MARS_BARYCENTER, "Mars barycenter.")
      .value("JUPITER_BARYCENTER", BodyId::JUPITER_BARYCENTER, "Jupiter barycenter.")
      .value("SATURN_BARYCENTER", BodyId::SATURN_BARYCENTER, "Saturn barycenter.")
      .value("URANUS_BARYCENTER", BodyId::URANUS_BARYCENTER, "Uranus barycenter.")
      .value("NEPTUNE_BARYCENTER", BodyId::NEPTUNE_BARYCENTER, "Neptune barycenter.")
      .value("PLUTO_BARYCENTER", BodyId::PLUTO_BARYCENTER, "Pluto barycenter.")
      .value("SUN", BodyId::SUN, "Sun.")
      .value("MERCURY", BodyId::MERCURY, "Mercury.")
      .value("VENUS", BodyId::VENUS, "Venus.")
      .value("EARTH", BodyId::EARTH, "Earth.")
      .value("MOON", BodyId::MOON, "Moon.")
      .value("MARS", BodyId::MARS, "Mars.")
      .value("PHOBOS", BodyId::PHOBOS, "Phobos.")
      .value("DEIMOS", BodyId::DEIMOS, "Deimos.")
      .value("JUPITER", BodyId::JUPITER, "Jupiter.")
      .export_values();

  py::enum_<Time>(m, "Time", "Time system / scale identifier.")
      .value("UT1", Time::UT1, "Universal Time 1.")
      .value("UTC", Time::UTC, "Coordinated Universal Time.")
      .value("TAI", Time::TAI, "International Atomic Time.")
      .value("TDB", Time::TDB, "Barycentric Dynamical Time.")
      .value("TT", Time::TT, "Terrestrial Time.")
      .value("TCG", Time::TCG, "Geocentric Coordinate Time.")
      .value("TCB", Time::TCB, "Barycentric Coordinate Time.")
      .value("GPS", Time::GPS, "GPS Time.")
      .value("TCL", Time::TCL, "Lunar Coordinate Time.")
      .value("LT", Time::LT, "Lunar Time.")
      .export_values();
}
