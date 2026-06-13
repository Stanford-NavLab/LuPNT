#pragma once

#include <string>
#include <tuple>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/graphs.h"
#include "lupnt/numerics/vector_macros.h"

namespace lupnt {

  extern const std::map<Time, std::string> time_to_string;
  extern const std::map<std::string, Time> string_to_time;

  Real ConvertTime(Real t, Time from, Time to);
  VecX ConvertTime(const VecX& t, Time from, Time to);
  bool IsCoordinateTimeScale(Time time);
  Real ConvertCoordinateTime(Real t, Time from, Time to);
  Real ConvertCoordinateTime(Real t, Time from, Time to, const Vec3& x_bcrs);
  VecX ConvertCoordinateTime(const VecX& t, Time from, Time to);
  VecX ConvertCoordinateTime(const VecX& t, Time from, Time to, const MatX3& x_bcrs);

  Real UtcToUt1(Real t_utc);
  Real Ut1ToUtc(Real t_ut1);

  Real TaiToUtc(Real t_tai);
  Real UtcToTai(Real t_utc);

  Real TaiToTt(Real t_tai);
  Real TtToTai(Real t_tt);

  Real TcgToTt(Real t_tcg);
  Real TtToTcg(Real t_tt);

  Real TtToTdb(Real t_tt);
  Real TtToTdb(Real t_tt, const Vec3& x_bcrs);
  Real TDBToTt(Real t_tdb);
  Real TDBToTt(Real t_tdb, const Vec3& x_bcrs);

  Real TaiToGps(Real t_tai);
  Real GpsToTai(Real t_gps);

  Real TcbToTdb(Real t_tcb);
  Real TtToTcb(Real t_tt);
  Real TtToTcb(Real t_tt, const Vec3& x_bcrs);
  Real TcbToTt(Real t_tcb, const Vec3& x_bcrs);

  Real MjdToTime(Real mjd);
  Real TimeToMjd(Real t);

  Real JdToTime(Real jd);
  Real TimeToJd(Real t);

  Real EarthRotationAngle(Real t_ut1);
  Real GregorianToMjd(int year, int month, int day, int hour = 0, int min = 0, Real sec = 0);
  Real GregorianToTime(int year, int month, int day, int hour = 0, int min = 0, Real sec = 0);
  Real GregorianToTime(const std::string& date);

  Real GreenwichMeanSiderealTime(Real mjd_ut1);
  Real GreenwichApparentSiderealTime(Real mjd_ut1);

  std::tuple<int, int, int, int, int, Real> MjdToGregorian(Real mjd);

  std::string MjdToGregorianString(Real mjd, int precision = 3);
  std::string TimeToGregorianString(Real t, int precision = 3);

  Real TcbToTcl(Real t_tcb);
  Real TcbToTcl(Real t_tcb, const Vec3& x_bcrs);
  VecX TcbToTcl(const VecX& t_tcb);
  VecX TcbToTcl(const VecX& t_tcb, const MatX3& x_bcrs);
  Real TclToTcb(Real t_tcl);
  Real TclToTcb(Real t_tcl, const Vec3& x_bcrs);
  VecX TclToTcb(const VecX& t_tcl);
  VecX TclToTcb(const VecX& t_tcl, const MatX3& x_bcrs);
  Real TclToLt(Real t_tcl);
  Real LtToTcl(Real t_lt);

  Real GetProperTimeCorrectionTcl(Real t_tcg, const Vec3& x_mci);
  VecX GetProperTimeCorrectionTcl(const VecX& t_tcg, const MatX& x_mci);

  // -----------------------------------------------------------------------
  // TL − TT conversions (Turyshev 2026, ApJ 997:97)
  // -----------------------------------------------------------------------

  /// TL − TT as a function of TDB epoch (Eq. 57).  Returns [s].
  Real TdbToLtMinusTt(Real t_tdb);

  /// TDB → TL direct conversion.
  Real TdbToLt(Real t_tdb);
  /// TL → TDB inversion (Newton-Raphson).
  Real LtToTdb(Real t_lt);

  /// TT → TL via TDB.
  Real TtToLt(Real t_tt);
  /// TL → TT inversion.
  Real LtToTt(Real t_lt);

  /// Fit TL−TT(TDB) over [t_start_tdb, t_end_tdb] with piecewise Chebyshev
  /// polynomials.  After this call TdbToLtMinusTt() is fast inside the window.
  void InitLtMinusTtFit(Real t_start_tdb, Real t_end_tdb,
                        double segment_length = 86400.0,  // 1 day [s]
                        int num_coeffs = 13);

  /// Clear the fitted model (revert to direct integration).
  void ClearLtMinusTtFit();

  /// Return true if a fitted model covers t_tdb.
  bool HasFittedLtMinusTt(Real t_tdb);

  VEC_DEF_REAL(UtcToUt1)
  VEC_DEF_REAL(Ut1ToUtc)
  VEC_DEF_REAL(TaiToUtc)
  VEC_DEF_REAL(UtcToTai)
  VEC_DEF_REAL(TaiToTt)
  VEC_DEF_REAL(TtToTai)
  VEC_DEF_REAL(TcgToTt)
  VEC_DEF_REAL(TtToTcg)
  VEC_DEF_REAL(TtToTdb)
  VEC_DEF_REAL(TDBToTt)
  VEC_DEF_REAL(TaiToGps)
  VEC_DEF_REAL(GpsToTai)
  VEC_DEF_REAL(TcbToTdb)
  VEC_DEF_REAL(TtToTcb)

  VEC_DEF_REAL(MjdToTime)
  VEC_DEF_REAL(TimeToMjd)
  VEC_DEF_REAL(JdToTime)
  VEC_DEF_REAL(TimeToJd)
  VEC_DEF_REAL(TcbToTcl)
  VEC_DEF_REAL(TclToTcb)
  VEC_DEF_REAL(TclToLt)
  VEC_DEF_REAL(LtToTcl)
  VEC_DEF_REAL(TdbToLtMinusTt)
  VEC_DEF_REAL(TdbToLt)
  VEC_DEF_REAL(LtToTdb)
  VEC_DEF_REAL(TtToLt)
  VEC_DEF_REAL(LtToTt)

  VEC_DEF_REAL(EarthRotationAngle)
  VEC_DEF_REAL(GreenwichMeanSiderealTime)
  VEC_DEF_REAL(GreenwichApparentSiderealTime)

}  // namespace lupnt
