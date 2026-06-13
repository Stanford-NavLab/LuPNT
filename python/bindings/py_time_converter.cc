// lupnt
#include <lupnt/conversions/coordinate_conversions.h>
#include <lupnt/conversions/time_conversions.h>
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

#include "lupnt/states/state.h"

// pybind11
#include <string>

#include "py_pybind11.h"
#include "py_vectorized_macros.h"

namespace py = pybind11;
using namespace lupnt;

void InitTimeConverter(py::module& m) {
  // ConvertTime
  m.def("convert_time", py::overload_cast<Real, Time, Time>(&ConvertTime), py::arg("t_tai"),
        py::arg("from_time"), py::arg("to_time"));
  m.def("convert_time", py::overload_cast<const VecX&, Time, Time>(&ConvertTime), py::arg("t_tai"),
        py::arg("from_time"), py::arg("to_time"));
  m.def("is_coordinate_time_scale", &IsCoordinateTimeScale, py::arg("time"));
  m.def("convert_coordinate_time", py::overload_cast<Real, Time, Time>(&ConvertCoordinateTime),
        py::arg("t"), py::arg("from_time"), py::arg("to_time"));
  m.def("convert_coordinate_time",
        py::overload_cast<Real, Time, Time, const Vec3&>(&ConvertCoordinateTime), py::arg("t"),
        py::arg("from_time"), py::arg("to_time"), py::arg("x_bcrs"));
  m.def("convert_coordinate_time",
        py::overload_cast<const VecX&, Time, Time>(&ConvertCoordinateTime), py::arg("t"),
        py::arg("from_time"), py::arg("to_time"));
  m.def("convert_coordinate_time",
        py::overload_cast<const VecX&, Time, Time, const MatX3&>(&ConvertCoordinateTime),
        py::arg("t"), py::arg("from_time"), py::arg("to_time"), py::arg("x_bcrs"));

  // GregorianToMjd
  m.def("gregorian_to_mjd", &GregorianToMjd, py::arg("year"), py::arg("month"), py::arg("day"),
        py::arg("hour"), py::arg("min"), py::arg("sec"));

  // GregorianToTime
  m.def("gregorian_to_time", py::overload_cast<int, int, int, int, int, Real>(&GregorianToTime),
        py::arg("year"), py::arg("month"), py::arg("day"), py::arg("hour"), py::arg("min"),
        py::arg("sec"));
  m.def("gregorian_to_time", py::overload_cast<const std::string&>(&GregorianToTime),
        py::arg("date"));

  // MjdToGregorian
  m.def("mjd_to_gregorian", py::overload_cast<Real>(&MjdToGregorian), py::arg("mjd"));

  m.def("mjd_to_gregorian_string", py::overload_cast<Real, int>(&MjdToGregorianString),
        py::arg("mjd"), py::arg("precision") = 3);
  m.def("time_to_gregorian_string", py::overload_cast<Real, int>(&TimeToGregorianString),
        py::arg("t"), py::arg("precision") = 3);

  // Lunar Related Conversions
  m.def("tcb_to_tcl", py::overload_cast<Real>(&TcbToTcl), py::arg("t_tcb"));
  m.def("tcb_to_tcl", py::overload_cast<Real, const Vec3&>(&TcbToTcl), py::arg("t_tcb"),
        py::arg("x_bcrs"));
  m.def("tcb_to_tcl", py::overload_cast<const VecX&>(&TcbToTcl), py::arg("t_tcb"));
  m.def("tcb_to_tcl", py::overload_cast<const VecX&, const MatX3&>(&TcbToTcl), py::arg("t_tcb"),
        py::arg("x_bcrs"));
  m.def("tcl_to_tcb", py::overload_cast<Real>(&TclToTcb), py::arg("t_tcl"));
  m.def("tcl_to_tcb", py::overload_cast<Real, const Vec3&>(&TclToTcb), py::arg("t_tcl"),
        py::arg("x_bcrs"));
  m.def("tcl_to_tcb", py::overload_cast<const VecX&>(&TclToTcb), py::arg("t_tcl"));
  m.def("tcl_to_tcb", py::overload_cast<const VecX&, const MatX3&>(&TclToTcb), py::arg("t_tcl"),
        py::arg("x_bcrs"));
  m.def("get_proper_time_correction_tcl",
        py::overload_cast<Real, const Vec3&>(&GetProperTimeCorrectionTcl), py::arg("t_tcg"),
        py::arg("x_mci"));
  m.def("get_proper_time_correction_tcl",
        py::overload_cast<const VecX&, const MatX&>(&GetProperTimeCorrectionTcl), py::arg("t_tcg"),
        py::arg("x_mci"));

  //   m.DEF_REAL("utc_to_ut1", UtcToUt1, "t_utc");
  //   m.DEF_REAL("ut1_to_utc", Ut1ToUtc, "t_ut1");
  //   m.DEF_REAL("tai_to_utc", TaiToUtc, "t_tai");
  //   m.DEF_REAL("utc_to_tai", UtcToTai, "t_utc");
  //   m.DEF_REAL("tai_to_tt", TaiToTt, "t_tai");
  //   m.DEF_REAL("tt_to_tai", TtToTai, "t_tt");
  //   m.DEF_REAL("tcg_to_tt", TcgToTt, "t_tcg");
  //   m.DEF_REAL("tt_to_tcg", TtToTcg, "t_tt");
  //   m.DEF_REAL("tt_to_tdb", TtToTdb, "t_tt");
  //   m.DEF_REAL("tdb_to_tt", TDBToTt, "t_tdb");
  //   m.DEF_REAL("tai_to_gps", TaiToGps, "t_tai");
  //   m.DEF_REAL("gps_to_tai", GpsToTai, "t_gps");
  //   m.DEF_REAL("tcb_to_tdb", TcbToTdb, "t_tcb");
  //   m.DEF_REAL("tt_to_tcb", TtToTcb, "t_tdb");

  m.DEF_REAL("mjd_to_time", MjdToTime, "mjd");
  m.DEF_REAL("time_to_mjd", TimeToMjd, "t");
  //   m.DEF_REAL("jd_to_time", JdToTime, "jd");
  //   m.DEF_REAL("time_to_jd", TimeToJd, "t");
}
