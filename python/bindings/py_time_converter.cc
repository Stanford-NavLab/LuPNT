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
        py::arg("from_time"), py::arg("to_time"),
        "Convert a time value from `from_time` to `to_time` (TAI, TT, TDB, UTC, UT1, GPS, ...) "
        "[s].");
  m.def("convert_time", py::overload_cast<const VecX&, Time, Time>(&ConvertTime), py::arg("t_tai"),
        py::arg("from_time"), py::arg("to_time"),
        "Vectorized: convert each time value from `from_time` to `to_time` [s].");
  m.def("is_coordinate_time_scale", &IsCoordinateTimeScale, py::arg("time"),
        "True if `time` is a coordinate time scale (TCB, TCG, TCL).");
  m.def("convert_coordinate_time", py::overload_cast<Real, Time, Time>(&ConvertCoordinateTime),
        py::arg("t"), py::arg("from_time"), py::arg("to_time"),
        "Convert between coordinate time scales (TCB/TCG/TCL) [s].");
  m.def("convert_coordinate_time",
        py::overload_cast<Real, Time, Time, const Vec3&>(&ConvertCoordinateTime), py::arg("t"),
        py::arg("from_time"), py::arg("to_time"), py::arg("x_bcrs"),
        "Convert between coordinate time scales using barycentric position `x_bcrs` for the "
        "position-dependent term [s].");
  m.def("convert_coordinate_time",
        py::overload_cast<const VecX&, Time, Time>(&ConvertCoordinateTime), py::arg("t"),
        py::arg("from_time"), py::arg("to_time"),
        "Vectorized: convert each value between coordinate time scales (TCB/TCG/TCL) [s].");
  m.def("convert_coordinate_time",
        py::overload_cast<const VecX&, Time, Time, const MatX3&>(&ConvertCoordinateTime),
        py::arg("t"), py::arg("from_time"), py::arg("to_time"), py::arg("x_bcrs"),
        "Vectorized: convert between coordinate time scales using per-epoch barycentric positions "
        "`x_bcrs` [s].");

  // GregorianToMjd
  m.def("gregorian_to_mjd", &GregorianToMjd, py::arg("year"), py::arg("month"), py::arg("day"),
        py::arg("hour"), py::arg("min"), py::arg("sec"),
        "Gregorian calendar date to Modified Julian Date [days].");

  // GregorianToTime
  m.def("gregorian_to_time", py::overload_cast<int, int, int, int, int, Real>(&GregorianToTime),
        py::arg("year"), py::arg("month"), py::arg("day"), py::arg("hour"), py::arg("min"),
        py::arg("sec"), "Gregorian calendar date to time [s past J2000].");
  m.def("gregorian_to_time", py::overload_cast<const std::string&>(&GregorianToTime),
        py::arg("date"), "Parse a calendar date string to time [s past J2000].");

  // MjdToGregorian
  m.def("mjd_to_gregorian", py::overload_cast<Real>(&MjdToGregorian), py::arg("mjd"),
        "Modified Julian Date to Gregorian (year, month, day, hour, min, sec).");

  m.def("mjd_to_gregorian_string", py::overload_cast<Real, int>(&MjdToGregorianString),
        py::arg("mjd"), py::arg("precision") = 3,
        "Format a Modified Julian Date as a Gregorian date string (`precision` fractional-second "
        "digits).");
  m.def("time_to_gregorian_string", py::overload_cast<Real, int>(&TimeToGregorianString),
        py::arg("t"), py::arg("precision") = 3,
        "Format a time [s] as a Gregorian date string (`precision` fractional-second digits).");

  // Lunar Related Conversions
  m.def("tcb_to_tcl", py::overload_cast<Real>(&TcbToTcl), py::arg("t_tcb"),
        "Barycentric Coordinate Time (TCB) to Lunar Coordinate Time (TCL) [s].");
  m.def("tcb_to_tcl", py::overload_cast<Real, const Vec3&>(&TcbToTcl), py::arg("t_tcb"),
        py::arg("x_bcrs"),
        "TCB to TCL using barycentric position `x_bcrs` for the position-dependent term [s].");
  m.def("tcb_to_tcl", py::overload_cast<const VecX&>(&TcbToTcl), py::arg("t_tcb"),
        "Vectorized: TCB to TCL [s].");
  m.def("tcb_to_tcl", py::overload_cast<const VecX&, const MatX3&>(&TcbToTcl), py::arg("t_tcb"),
        py::arg("x_bcrs"), "Vectorized: TCB to TCL using per-epoch barycentric positions [s].");
  m.def("tcl_to_tcb", py::overload_cast<Real>(&TclToTcb), py::arg("t_tcl"),
        "Lunar Coordinate Time (TCL) to Barycentric Coordinate Time (TCB) [s].");
  m.def("tcl_to_tcb", py::overload_cast<Real, const Vec3&>(&TclToTcb), py::arg("t_tcl"),
        py::arg("x_bcrs"),
        "TCL to TCB using barycentric position `x_bcrs` for the position-dependent term [s].");
  m.def("tcl_to_tcb", py::overload_cast<const VecX&>(&TclToTcb), py::arg("t_tcl"),
        "Vectorized: TCL to TCB [s].");
  m.def("tcl_to_tcb", py::overload_cast<const VecX&, const MatX3&>(&TclToTcb), py::arg("t_tcl"),
        py::arg("x_bcrs"), "Vectorized: TCL to TCB using per-epoch barycentric positions [s].");
  m.def("get_proper_time_correction_tcl",
        py::overload_cast<Real, const Vec3&>(&GetProperTimeCorrectionTcl), py::arg("t_tcg"),
        py::arg("x_mci"),
        "Proper-time correction for Lunar Coordinate Time (TCL) at Moon-centered inertial position "
        "`x_mci` [s].");
  m.def("get_proper_time_correction_tcl",
        py::overload_cast<const VecX&, const MatX&>(&GetProperTimeCorrectionTcl), py::arg("t_tcg"),
        py::arg("x_mci"),
        "Vectorized: TCL proper-time correction at per-epoch Moon-centered inertial positions "
        "`x_mci` [s].");

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

  // ---------------------------------------------------------------------
  // TT <-> TDB models
  // ---------------------------------------------------------------------
  py::enum_<TtTdbModel>(m, "TtTdbModel", "TT<->TDB computation selected by set_tt_tdb_model.")
      .value("ANALYTIC", TtTdbModel::ANALYTIC,
             "Two-term analytic series (fast, ~30 us vs DE440t). Always returns the series.")
      .value("DE440_INTEGRAL", TtTdbModel::DE440_INTEGRAL,
             "DE440 Eq. (3) relativistic integral (Park et al. 2021); integrates from T_0.")
      .value("FITTED", TtTdbModel::FITTED,
             "Piecewise Chebyshev fit of DE440t (default). Uses an initialised or "
             "auto-fit model; falls back to the series only if none can be built.");

  m.def("set_tt_tdb_model", &SetTtTdbModel, py::arg("model"),
        "Select the TT<->TDB computation (ANALYTIC / DE440_INTEGRAL / FITTED).");
  m.def("get_tt_tdb_model", &GetTtTdbModel, "Current TT<->TDB model.");
  m.def("set_de440_tt_tdb_step", &SetDe440TtTdbStep, py::arg("step_s"),
        "Trapezoidal step [s] for the DE440 Eq. (3) integral (default 864 s = 0.01 day).");
  m.def("get_de440_tt_tdb_step", &GetDe440TtTdbStep,
        "Trapezoidal step [s] for the DE440 Eq. (3) integral.");

  m.def("set_tt_tdb_auto_fit", &SetTtTdbAutoFit, py::arg("enable"),
        "Build a DE440t Chebyshev fit on demand when none covers the requested epoch "
        "(default: True). Makes the DEFAULT TT<->TDB accuracy ~0.4 ps instead of the "
        "~17 us analytic series, for convert_time, Epoch and the offset APIs alike. "
        "Set False for the historical behaviour or to avoid touching SPICE.");
  m.def("get_tt_tdb_auto_fit", &GetTtTdbAutoFit, "Whether TT<->TDB auto-fitting is enabled.");

  m.def("init_tdb_minus_tcl_fit", &InitTdbMinusTclFit, py::arg("t_start_tdb"), py::arg("t_end_tdb"),
        py::arg("segment_length") = 4.0 * 86400.0, py::arg("num_coeffs") = 13,
        "Fit TDB-TCL over a window with piecewise Chebyshev polynomials.");
  m.def("clear_tdb_minus_tcl_fit", &ClearTdbMinusTclFit, "Clear the fitted TDB-TCL model.");
  m.def("has_fitted_tdb_minus_tcl", &HasFittedTdbMinusTcl, py::arg("t_tdb"),
        "Whether a fitted TDB-TCL model covers this epoch.");
  m.def("set_tdb_tcl_auto_fit", &SetTdbTclAutoFit, py::arg("enable"),
        "Build a TDB-TCL Chebyshev fit on demand (default: True). Without it, each "
        "TdbMinusTcl call integrates from T_0 (1977), ~12 s per call.");
  m.def("get_tdb_tcl_auto_fit", &GetTdbTclAutoFit, "Whether TDB<->TCL auto-fitting is enabled.");

  m.def("tdb_minus_tt_de440", py::overload_cast<Real>(&TdbMinusTtDe440), py::arg("t_tdb"),
        "TDB - TT [s] at the geocenter from DE440 Eq. (3) (Park et al. 2021, AJ 161:105).");
  m.def("tdb_minus_tt_de440", py::overload_cast<const VecX&>(&TdbMinusTtDe440), py::arg("t_tdb"),
        "Vectorized geocentric TDB - TT [s] from DE440 Eq. (3). Uses one sorted sweep of the "
        "integral instead of re-integrating from T_0 per epoch -- strongly preferred for arrays.");
  m.def("tdb_minus_tt_de440", py::overload_cast<Real, const Vec3&>(&TdbMinusTtDe440),
        py::arg("t_tdb"), py::arg("x_bcrs"),
        "TDB - TT [s] from DE440 Eq. (3) for a station at BCRS position `x_bcrs` [m].");

  // --- Offset accessors: avoid the absolute-epoch ULP floor ----------------
  // At |t| ~ 1e9 s one ULP is ~0.25 us, so `convert(t) - t` snaps a small
  // time-scale difference to that grid. These return the difference directly.
  m.def("tt_minus_tdb", py::overload_cast<Real>(&TtMinusTdb), py::arg("t_tdb"),
        "TT - TDB [s] under the active model/fit, as a full-precision offset "
        "(no absolute epoch is formed, so no ~0.25 us ULP floor).");
  m.def("tt_minus_tdb", py::overload_cast<const VecX&>(&TtMinusTdb), py::arg("t_tdb"),
        "Vectorized TT - TDB [s] offset.");
  m.def("tdb_minus_tcl", py::overload_cast<Real>(&TdbMinusTcl), py::arg("t_tdb"),
        "TDB - TCL [s] at the lunar centre, as a full-precision offset.");
  m.def("tdb_minus_tcl", py::overload_cast<const VecX&>(&TdbMinusTcl), py::arg("t_tdb"),
        "Vectorized TDB - TCL [s] offset.");
  m.def("tdb_minus_lt", py::overload_cast<Real>(&TdbMinusLt), py::arg("t_tdb"),
        "TDB - TL [s] as a full-precision offset (uses the active TT<->TDB model).");
  m.def("tdb_minus_lt", py::overload_cast<const VecX&>(&TdbMinusLt), py::arg("t_tdb"),
        "Vectorized TDB - TL [s] offset.");

  m.def("tdb_to_tt", py::overload_cast<Real>(&TDBToTt), py::arg("t_tdb"),
        "Convert TDB to TT [s] using the active TT<->TDB model / Chebyshev fit.");
  m.def("tt_to_tdb", py::overload_cast<Real>(&TtToTdb), py::arg("t_tt"),
        "Convert TT to TDB [s] using the active TT<->TDB model / Chebyshev fit.");

  // ---------------------------------------------------------------------
  // TT - TDB Chebyshev fit (sampled from the DE440t TT-TDB ephemeris)
  // ---------------------------------------------------------------------
  m.def("init_tt_minus_tdb_fit", &InitTtMinusTdbFit, py::arg("t_start_tdb"), py::arg("t_end_tdb"),
        py::arg("segment_length") = 16.0 * 86400.0, py::arg("num_coeffs") = 13,
        "Fit TT-TDB(TDB) over a window with piecewise Chebyshev polynomials sampled from the "
        "DE440t TT-TDB ephemeris. Afterwards tdb_to_tt/tt_to_tdb evaluate the fit (no SPICE).");
  m.def("clear_tt_minus_tdb_fit", &ClearTtMinusTdbFit,
        "Clear the fitted TT-TDB model (revert to the analytic series).");
  m.def("has_fitted_tt_minus_tdb", &HasFittedTtMinusTdb, py::arg("t"),
        "True if a fitted TT-TDB model covers epoch `t`.");

  // TL - TT fit controls (Turyshev 2026 Eq. 57)
  m.def("init_lt_minus_tt_fit", &InitLtMinusTtFit, py::arg("t_start_tdb"), py::arg("t_end_tdb"),
        py::arg("segment_length") = 86400.0, py::arg("num_coeffs") = 13,
        "Fit TL-TT(TDB) over a window with piecewise Chebyshev polynomials.");
  m.def("clear_lt_minus_tt_fit", &ClearLtMinusTtFit, "Clear the fitted TL-TT model.");
  m.def("tdb_to_lt_minus_tt", py::overload_cast<Real>(&TdbToLtMinusTt), py::arg("t_tdb"),
        "TL - TT [s] as a function of TDB (Turyshev 2026 Eq. 57).");

  m.DEF_REAL("mjd_to_time", MjdToTime, "mjd");
  m.DEF_REAL("time_to_mjd", TimeToMjd, "t");
  //   m.DEF_REAL("jd_to_time", JdToTime, "jd");
  //   m.DEF_REAL("time_to_jd", TimeToJd, "t");
}
