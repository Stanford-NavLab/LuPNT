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

  // -----------------------------------------------------------------------
  // TT <-> TDB model selection
  // -----------------------------------------------------------------------

  /// Which model TDBToTt()/TtToTdb() use when no Chebyshev fit (see
  /// InitTtMinusTdbFit) covers the requested epoch.
  enum class TtTdbModel {
    /// Two-term analytic series. Fast; differs from the integrated DE440t
    /// TT-TDB ephemeris by up to ~30 us. Selecting this always returns the
    /// series, ignoring any fit.
    ANALYTIC,
    /// DE440 Eq. (3) relativistic integral (Park et al. 2021, AJ 161:105).
    /// Much slower -- integrates from T_0 on every call -- but reproduces the
    /// JPL TT-TDB relation directly from the ephemeris.
    DE440_INTEGRAL,
    /// Piecewise Chebyshev fit of the DE440t ephemeris (default). Uses a fit
    /// covering the epoch if one has been initialised (InitTtMinusTdbFit), else
    /// builds one on demand when auto-fit is enabled (SetTtTdbAutoFit, on by
    /// default); falls back to the series only if no fit can be produced.
    /// Sub-nanosecond against DE440t.
    FITTED,
  };

  void SetTtTdbModel(TtTdbModel model);
  TtTdbModel GetTtTdbModel();

  /// Build a DE440t Chebyshev fit on demand when no fit covers the requested
  /// epoch (default: ON). This makes the DEFAULT TT<->TDB accuracy ~0.4 ps
  /// instead of the ~17 us of the legacy analytic series, for every caller --
  /// ConvertTime, Epoch, and the offset APIs alike.
  ///
  /// The fit is built once per decade-wide window (snapped to a fixed grid, so
  /// sequential epochs reuse it) and costs a few thousand SPICE reads. If the
  /// DE440t TT-TDB segment is unavailable the failure is recorded once and all
  /// callers fall back to the analytic series.
  ///
  /// PERFORMANCE: only ONE window is cached at a time. Sequential or clustered
  /// access is cheap, but repeatedly jumping between epochs in DIFFERENT decade
  /// windows rebuilds the fit on every crossing (~15 s each). For a workload
  /// that samples a multi-decade span non-monotonically, call InitTtMinusTdbFit
  /// (and InitTdbMinusTclFit) once over the whole span up front instead.
  ///
  /// Set false for the historical behaviour, or for runs that must not touch
  /// SPICE.
  void SetTtTdbAutoFit(bool enable);
  bool GetTtTdbAutoFit();

  /// Trapezoidal step [s] used by the DE440 Eq. (3) integral (default 0.01 day).
  void SetDe440TtTdbStep(double step_s);
  double GetDe440TtTdbStep();

  // -----------------------------------------------------------------------
  // Offset accessors (avoid the absolute-epoch ULP floor)
  //
  // LuPNT epochs are doubles holding seconds from J2000, so at present-day
  // dates (|t| ~ 1e9 s) one ULP is ~0.25 us. Computing a small time-scale
  // difference as `Convert(t) - t` therefore snaps it to that grid. These
  // accessors return the difference directly, never forming an absolute
  // epoch, and so retain full double precision (sub-femtosecond).
  // -----------------------------------------------------------------------

  /// TT - TDB [s] under the active TT<->TDB model / Chebyshev fit.
  Real TtMinusTdb(Real t_tdb);
  /// TDB - TCL [s] at the lunar centre.
  Real TdbMinusTcl(Real t_tdb);
  /// TDB - TL [s] (uses the active TT<->TDB model).
  Real TdbMinusLt(Real t_tdb);

  /// TDB - TT [s] from DE440 Eq. (3) at the geocenter.
  ///
  /// RESIDUAL MODEL ERROR: ~0.14 ns/yr secular drift relative to the DE440t
  /// TT-TDB ephemeris (measured 2020-2035); ~7 ns rms. The periodic structure is
  /// reproduced to 0.006 ns rms -- only the rate is imperfect.
  ///
  /// w_0E sums the 10 ephemeris bodies plus two small-body populations that DE440
  /// integrates discretely but for which LuPNT has no ephemerides: the main
  /// asteroid belt (343 bodies in DE440) and the Kuiper belt (30 discrete KBOs
  /// plus a 36-point ring at 44 au). Each is modelled as a uniform circular ring
  /// via RingPotential(), using PUBLISHED masses -- see GM_ASTEROID_BELT and
  /// GM_KUIPER_BELT in constants.h. Their combined potential accounts for ~80%
  /// of the secular difference against DE440t, consistent with the ~77% those
  /// published masses predict.
  ///
  /// The remaining ~0.14 ns/yr is ~0.40 m^2/s^2 of potential, within ~2 sigma of
  /// the published mass uncertainties (Kuiper belt alone is +-15%), the
  /// single-radius ring approximation, and the fact that DE440 fitted its own
  /// ring mass rather than adopting the published total.
  ///
  /// Do NOT tune GM_ASTEROID_BELT / GM_KUIPER_BELT to drive this to zero -- they
  /// are literature values and tuning them reproduces the number without the
  /// physics. For sub-ns absolute agreement with JPL use InitTtMinusTdbFit()
  /// (0.0004 ns vs DE440t).
  ///
  /// Ruled out as contributors, by measurement:
  ///   * quadrature -- refining 864s -> 216s moves the drift by 0.0005 ns/yr
  ///   * the c^-4 terms -- zeroing them shifts it by exactly the predicted
  ///     +3.46 ns/yr those terms supply
  ///   * body IDs / GMs -- all resolve to the correct barycentres with system GMs
  ///   * L_B / L_G / L_L / L_H / L_M -- audited against Turyshev et al. 2025 Table 2
  Real TdbMinusTtDe440(Real t_tdb);
  /// Vectorized geocentric DE440 Eq. (3) TDB - TT [s]. Uses a single sorted
  /// sweep of the integral instead of re-integrating from T_0 per epoch.
  VecX TdbMinusTtDe440(const VecX& t_tdb);
  /// TDB - TT [s] from DE440 Eq. (3) for a station at BCRS position `x_bcrs`.
  Real TdbMinusTtDe440(Real t_tdb, const Vec3& x_bcrs);

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

  // -----------------------------------------------------------------------
  // TT − TDB high-fidelity fit (DE440t TT-TDB ephemeris)
  // -----------------------------------------------------------------------

  /// Fit TT−TDB(TDB) over [t_start_tdb, t_end_tdb] with piecewise Chebyshev
  /// polynomials sampled from the DE440t TT-TDB ephemeris (via
  /// spice::ConvertTime).  After this call TDBToTt()/TtToTdb() use the fit
  /// (no SPICE call) inside the window; outside it they fall back to the
  /// analytic series.
  void InitTtMinusTdbFit(Real t_start_tdb, Real t_end_tdb,
                         double segment_length = 16.0 * 86400.0,  // 16 days [s]
                         int num_coeffs = 13);

  /// Clear the fitted TT−TDB model (revert to the analytic series).
  void ClearTtMinusTdbFit();

  // -----------------------------------------------------------------------
  // TDB − TCL fit
  //
  // TdbMinusTcl() integrates from T_0 (1977) on every call (~12 s), so a fit is
  // required for any per-step lunar time conversion. Auto-fitting is ON by
  // default and builds one decade-wide window on demand.
  // -----------------------------------------------------------------------

  void InitTdbMinusTclFit(Real t_start_tdb, Real t_end_tdb, double segment_length = 4.0 * 86400.0,
                          int num_coeffs = 13);
  void ClearTdbMinusTclFit();
  bool HasFittedTdbMinusTcl(Real t_tdb);
  void SetTdbTclAutoFit(bool enable);
  bool GetTdbTclAutoFit();

  /// Return true if a fitted TT−TDB model covers epoch t.
  bool HasFittedTtMinusTdb(Real t);

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
  VEC_DEF_REAL(TtMinusTdb)
  VEC_DEF_REAL(TdbMinusTcl)
  VEC_DEF_REAL(TdbMinusLt)
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
