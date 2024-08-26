// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coordinates.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/time_converter.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <string>

#include "py_vectorized_macros.cc"

namespace py = pybind11;
using namespace lupnt;

void init_time_converter(py::module &m) {
  // ConvertTime
  m.def(
      "convert_time",
      [](double t, const Time &from, const Time &to) -> double {
        return ConvertTime(t, from, to).val();
      },
      py::arg("t"), py::arg("from"), py::arg("to"));
  m.def(
      "convert_time",
      [](const VecXd &t, const Time &from, const Time &to) -> VecXd {
        return ConvertTime(t.cast<Real>().eval(), from, to).cast<double>();
      },
      py::arg("t"), py::arg("from"), py::arg("to"));

  // Gregorian2MJD
  m.def(
      "gregorian2mjd",
      [](int year, int month, int day, int hour = 0, int min = 0, double sec = 0) -> double {
        return Gregorian2MJD(year, month, day, hour, min, sec).val();
      },
      py::arg("year"), py::arg("month"), py::arg("day"), py::arg("hour") = 0, py::arg("min") = 0,
      py::arg("sec") = 0);

  // Gregorian2Time
  m.def(
      "gregorian2time",
      [](int year, int month, int day, int hour = 0, int min = 0, double sec = 0) -> double {
        return Gregorian2Time(year, month, day, hour, min, sec).val();
      },
      py::arg("year"), py::arg("month"), py::arg("day"), py::arg("hour") = 0, py::arg("min") = 0,
      py::arg("sec") = 0);

  // MJD2Gregorian
  m.def(
      "mjd2gregorian",
      [](double mjd) -> std::tuple<int, int, int, int, int, double> {
        auto [year, month, day, hour, min, sec] = MJD2Gregorian(mjd);
        return std::make_tuple(year, month, day, hour, min, sec.val());
      },
      py::arg("mjd"));

  m.def(
      "mjd2gregorian_string",
      [](double mjd, int precision) -> std::string { return MJD2GregorianString(mjd, precision); },
      py::arg("mjd"), py::arg("precision") = 3);
  m.def(
      "time2gregorian_string",
      [](double t, int precision) -> std::string { return Time2GregorianString(t, precision); },
      py::arg("t"), py::arg("precision") = 3);

  VEC_BIND_REAL("utc2ut1", UTC2UT1, "t_utc")
  VEC_BIND_REAL("ut12utc", UT12UTC, "t_ut1")
  VEC_BIND_REAL("tai2utc", TAI2UTC, "t_tai")
  VEC_BIND_REAL("utc2tai", UTC2TAI, "t_utc")
  VEC_BIND_REAL("tai2tt", TAI2TT, "t_tai")
  VEC_BIND_REAL("tt2tai", TT2TAI, "t_tt")
  VEC_BIND_REAL("tcg2tt", TCG2TT, "t_tcg")
  VEC_BIND_REAL("tt2tcg", TT2TCG, "t_tt")
  VEC_BIND_REAL("tt2tdb", TT2TDB, "t_tt")
  VEC_BIND_REAL("tdb2tt", TDB2TT, "t_tdb")
  VEC_BIND_REAL("tai2gps", TAI2GPS, "t_tai")
  VEC_BIND_REAL("gps2tai", GPS2TAI, "t_gps")
  VEC_BIND_REAL("tcb2tdb", TCB2TDB, "t_tcb")
  VEC_BIND_REAL("tt2tcb", TT2TCB, "t_tdb")
  VEC_BIND_REAL("mjd2time", MJD2Time, "mjd")
  VEC_BIND_REAL("time2mjd", Time2MJD, "t")
  VEC_BIND_REAL("jd2time", JD2Time, "jd")
  VEC_BIND_REAL("time2jd", Time2JD, "t")
}
