#include "lupnt/conversions/epoch.h"

#include <cmath>
#include <cstdio>

#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/error.h"
#include "lupnt/interfaces/eop.h"
#include "lupnt/interfaces/tai_utc.h"
#include "lupnt/numerics/graphs.h"

namespace lupnt {

  void Epoch::Normalize() {
    // floor() of the value part; subtracting an integer keeps the derivative.
    double f = frac_.val();
    if (!std::isfinite(f)) return;
    double k = std::floor(f);
    if (k != 0.0) {
      sec_ += static_cast<int64_t>(k);
      frac_ -= k;
    }
  }

  Epoch Epoch::FromSeconds(Real t, Time scale) {
    double v = t.val();
    double k = std::floor(v);
    // Keep the fractional part as a Real so any derivative on `t` survives.
    return Epoch(static_cast<int64_t>(k), t - k, scale);
  }

  namespace {
    /// Days from 1970-01-01 in the proleptic Gregorian calendar, exactly.
    /// Howard Hinnant's civil-from-days algorithm; integer arithmetic only.
    int64_t DaysFromCivil(int64_t y, int64_t m, int64_t d) {
      y -= (m <= 2);
      const int64_t era = (y >= 0 ? y : y - 399) / 400;
      const int64_t yoe = y - era * 400;  // [0, 399]
      const int64_t doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
      const int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;  // [0, 146096]
      return era * 146097 + doe - 719468;
    }

    // 1970-01-01 -> 2000-01-01 is 10957 days; the J2000 origin is at 12:00, so a
    // civil midnight sits 43200 s before the day's J2000 second count.
    constexpr int64_t kDaysCivilToJ2000 = 10957;
    constexpr int64_t kJ2000NoonOffset = 43200;

    // Before the Gregorian reform GregorianToMjd switches to the Julian calendar;
    // keep that behaviour rather than silently reinterpreting historical dates.
    bool IsGregorianEra(int year, int month, int day) {
      return (10000L * year + 100L * month + day) > 15821004L;
    }
  }  // namespace

  /// @brief Build an Epoch from a calendar date, exactly.
  ///
  /// The integer part is computed with integer arithmetic and is therefore
  /// exact; only the sub-second remainder is carried as a Real.
  ///
  /// Routing a calendar date through a Modified Julian Date instead would round
  /// it: MJD counts days from 1858, so at present-day epochs one ULP of an MJD
  /// is ~1.1 us -- 8x coarser than seconds-from-J2000, and the rounding lands on
  /// the largest quantity in the chain. Whole minutes survive that (their
  /// fraction-of-day is a dyadic rational) but arbitrary clock times do not:
  /// 2020-07-04T12:34:56 comes out 238 ns late, 2020-01-01T12:34:56.789 comes
  /// out 358 ns early.
  Epoch Epoch::FromGregorian(int year, int month, int day, int hour, int min, Real sec,
                             Time scale) {
    if (!IsGregorianEra(year, month, day)) {
      // Julian-calendar dates: defer to the MJD path, which handles the reform.
      double s_int = std::floor(sec.val());
      Epoch e = FromSeconds(GregorianToTime(year, month, day, hour, min, Real(s_int)), scale);
      return e + (sec - s_int);
    }
    const int64_t days = DaysFromCivil(year, month, day) - kDaysCivilToJ2000;
    const double s = sec.val();
    const double s_floor = std::floor(s);
    const int64_t whole = days * 86400 + static_cast<int64_t>(hour) * 3600
                          + static_cast<int64_t>(min) * 60 + static_cast<int64_t>(s_floor)
                          - kJ2000NoonOffset;
    return Epoch(whole, sec - s_floor, scale);
  }

  Epoch Epoch::FromGregorian(const std::string& date, Time scale) {
    int year, month, day, hour, min;
    double sec;
    if (std::sscanf(date.c_str(), "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &min, &sec)
        == 6) {
      return FromGregorian(year, month, day, hour, min, Real(sec), scale);
    }
    return FromSeconds(GregorianToTime(date), scale);
  }

  Real Epoch::ToMjd() const { return TimeToMjd(ToSeconds()); }

  std::string Epoch::ToGregorianString(int precision) const {
    return TimeToGregorianString(ToSeconds(), precision);
  }

  Real Epoch::operator-(const Epoch& other) const {
    LUPNT_CHECK(scale_ == other.scale_,
                fmt::format("Cannot subtract epochs in different time scales ({} and {}); "
                            "convert one with Epoch::To() first",
                            time_to_string.at(scale_), time_to_string.at(other.scale_)),
                "Epoch");
    // The integer difference is exact in int64 and converts exactly to double
    // for any |difference| < 2^53 s (~285 million years).
    double whole = static_cast<double>(sec_ - other.sec_);
    return Real(whole) + (frac_ - other.frac_);
  }

  bool Epoch::operator==(const Epoch& o) const {
    return scale_ == o.scale_ && sec_ == o.sec_ && frac_.val() == o.frac_.val();
  }

  bool Epoch::operator<(const Epoch& o) const {
    LUPNT_CHECK(scale_ == o.scale_,
                fmt::format("Cannot order epochs in different time scales ({} and {})",
                            time_to_string.at(scale_), time_to_string.at(o.scale_)),
                "Epoch");
    if (sec_ != o.sec_) return sec_ < o.sec_;
    return frac_.val() < o.frac_.val();
  }

  namespace {

    /// Offset `(scale - TAI)` [s] at the instant given by `t`, **which must be
    /// a TAI reading** (seconds from J2000).
    ///
    /// Every supported scale sits within ~100 s of TAI over the modern era
    /// (leap seconds ~37 s, TT-TAI = 32.184 s, TCB-TT ~ 21 s by 2035), so all
    /// of these are exact to ~1e-14 s in a double. Building conversions out of
    /// these small offsets -- rather than differencing absolute epochs -- is
    /// what keeps `Epoch::To` sub-nanosecond.
    ///
    /// The argument must be TAI, not "whatever scale the epoch happens to be
    /// in": passing a TDB reading here instead would shift the evaluation
    /// instant by ~32 s, which through d(TT-TDB)/dt ~ 3.2e-10 s/s injects a
    /// ~10 ns error -- above the precision this class exists to provide.
    /// `TimeScaleOffset` is responsible for establishing that TAI reference.
    ///
    /// Using the (lossy) absolute seconds for `t` is harmless: the offsets
    /// vary by <1e-8 s per second, so the ~2.5e-7 s ULP of the argument
    /// perturbs the result by <1e-15 s.
    // ------------------------------------------------------------------
    // Time-scale conversion graph
    // ------------------------------------------------------------------
    //
    // Each edge returns the OFFSET (to - from) [s] given the epoch's reading in
    // the `from` scale. A conversion walks the shortest registered route and
    // sums those offsets, so two large absolute epochs are never differenced and
    // the result keeps full double precision rather than the ~245 ns epoch ULP.
    //
    // Edge costs span orders of magnitude:
    //
    //     constant   TAI<->TT, TAI<->GPS
    //     table      TAI<->UTC (leap seconds), UTC<->UT1 (EOP)
    //     linear     TT<->TCG, TDB<->TCB, TCL<->LT   (a rescaling by L ~ 1e-8)
    //     model      TT<->TDB                        (Chebyshev fit / series)
    //     integral   TDB<->TCL                       (integrates from T_0, 1977)
    //
    // Routing matters because of that spread: TCL <-> LT and TDB <-> TCB are
    // single linear edges and must not be reached via a path that drags in the
    // TDB <-> TCL integral.
    //
    // The search minimises hops, which coincides with minimum cost for this edge
    // set: the one expensive edge (TDB <-> TCL) is also the only bridge to the
    // lunar scales, so no cheaper detour around it exists. Adding an edge that
    // breaks that property would require a cost-weighted search.
    using TimeOffsetFn = std::function<Real(Real)>;

    Real T0Seconds() { return MjdToTime(MJD_COORDINATE_TT_TCG_TCB); }

    const std::map<std::pair<Time, Time>, TimeOffsetFn>& TimeEdges() {
      static const std::map<std::pair<Time, Time>, TimeOffsetFn> edges = {
          // --- constant offsets ---
          {{Time::TAI, Time::TT}, [](Real) { return Real(TT_TAI_OFFSET); }},
          {{Time::TT, Time::TAI}, [](Real) { return Real(-TT_TAI_OFFSET); }},
          {{Time::TAI, Time::GPS}, [](Real) { return Real(-19.0); }},
          {{Time::GPS, Time::TAI}, [](Real) { return Real(19.0); }},

          // --- leap seconds / Earth orientation ---
          {{Time::TAI, Time::UTC},
           [](Real t_tai) { return Real(-GetTaiUtcDifference(TimeToMjd(t_tai).val())); }},
          {{Time::UTC, Time::TAI},
           [](Real t_utc) { return Real(GetTaiUtcDifference(TimeToMjd(t_utc).val())); }},
          {{Time::UTC, Time::UT1},
           [](Real t_utc) { return GetUt1UtcDifference(TimeToMjd(t_utc)); }},
          {{Time::UT1, Time::UTC},
           [](Real t_ut1) { return -GetUt1UtcDifference(TimeToMjd(t_ut1)); }},

          // --- linear rescalings (large elapsed time only ever scaled by L ~ 1e-8) ---
          {{Time::TT, Time::TCG},
           [](Real t_tt) { return L_G / (1.0 - L_G) * (t_tt - T0Seconds()); }},
          {{Time::TCG, Time::TT}, [](Real t_tcg) { return -L_G * (t_tcg - T0Seconds()); }},
          {{Time::TDB, Time::TCB},
           [](Real t_tdb) { return (L_B * (t_tdb - T0Seconds()) - TDB_0) / (1.0 - L_B); }},
          {{Time::TCB, Time::TDB}, [](Real t_tcb) { return -L_B * (t_tcb - T0Seconds()) + TDB_0; }},
          {{Time::TCL, Time::LT}, [](Real t_tcl) { return -L_L * (t_tcl - T0Seconds()); }},
          {{Time::LT, Time::TCL},
           [](Real t_lt) { return L_L / (1.0 - L_L) * (t_lt - T0Seconds()); }},

          // --- TT <-> TDB: the active model (Chebyshev fit / DE440 integral / series) ---
          {{Time::TT, Time::TDB}, [](Real t_tt) { return -TtMinusTdb(t_tt); }},
          {{Time::TDB, Time::TT}, [](Real t_tdb) { return TtMinusTdb(t_tdb); }},

          // --- TDB <-> TCL: the only genuinely expensive edge (integrates from T_0) ---
          {{Time::TDB, Time::TCL}, [](Real t_tdb) { return -TdbMinusTcl(t_tdb); }},
          {{Time::TCL, Time::TDB}, [](Real t_tcl) { return TdbMinusTcl(t_tcl); }},
      };
      return edges;
    }

  }  // namespace

  Real TimeScaleOffset(const Epoch& e, Time to) {
    if (e.scale() == to) return 0.0;

    const auto& edges = TimeEdges();
    LUPNT_CHECK(
        std::any_of(edges.begin(), edges.end(),
                    [&](const auto& kv) { return kv.first.first == e.scale(); }),
        fmt::format("Epoch does not support the time scale '{}'", time_to_string.at(e.scale())),
        "Epoch");
    LUPNT_CHECK(std::any_of(edges.begin(), edges.end(),
                            [&](const auto& kv) { return kv.first.second == to; }),
                fmt::format("Epoch does not support the time scale '{}'", time_to_string.at(to)),
                "Epoch");

    // Walk the cheapest registered route, accumulating OFFSETS. `total` only
    // ever sums small quantities, so the result keeps full double precision;
    // `t` tracks the reading in the current scale purely as the argument for the
    // next edge, where a ~245 ns error costs <1e-15 s because every edge varies
    // by <1e-8 s per second.
    const std::vector<Time> path = FindShortestPath(e.scale(), to, edges);
    Real t = e.ToSeconds();
    Real total = 0.0;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
      const Real off = edges.at({path[i], path[i + 1]})(t);
      total += off;
      t += off;
    }
    return total;
  }

  Epoch Epoch::To(Time target) const {
    if (target == scale_) return *this;
    Real offset = TimeScaleOffset(*this, target);
    return Epoch(sec_, frac_ + offset, target);
  }

  // =========================================================================
  // EpochSeries
  // =========================================================================

  void EpochSeries::Normalize() {
    for (int i = 0; i < size(); i++) {
      double f = frac_(i).val();
      if (!std::isfinite(f)) continue;
      double k = std::floor(f);
      if (k != 0.0) {
        sec_[i] += static_cast<int64_t>(k);
        frac_(i) -= k;
      }
    }
  }

  EpochSeries EpochSeries::Linspace(const Epoch& start, Real step, int n) {
    LUPNT_CHECK(n >= 0, "EpochSeries::Linspace: n must be non-negative", "EpochSeries");
    std::vector<int64_t> sec(n);
    VecX frac(n);

    // Split the step into an exact integer part and a small remainder, then
    // accumulate the integer part in int64. This keeps a long grid exact:
    // `start + i*step` never goes through a rounded large double.
    const double step_v = step.val();
    const double step_int = std::floor(step_v);
    const Real step_frac = step - step_int;

    for (int i = 0; i < n; i++) {
      sec[i] = start.seconds() + static_cast<int64_t>(step_int) * i;
      frac(i) = start.fraction() + step_frac * i;
    }
    return EpochSeries(std::move(sec), std::move(frac), start.scale());
  }

  EpochSeries EpochSeries::FromSeconds(const VecX& t, Time scale) {
    const int n = static_cast<int>(t.size());
    std::vector<int64_t> sec(n);
    VecX frac(n);
    for (int i = 0; i < n; i++) {
      double k = std::floor(t(i).val());
      sec[i] = static_cast<int64_t>(k);
      frac(i) = t(i) - k;
    }
    return EpochSeries(std::move(sec), std::move(frac), scale);
  }

  VecX EpochSeries::ToSeconds() const {
    VecX out(size());
    for (int i = 0; i < size(); i++) out(i) = Real(static_cast<double>(sec_[i])) + frac_(i);
    return out;
  }

  VecX EpochSeries::Since(const Epoch& ref) const {
    LUPNT_CHECK(ref.scale() == scale_,
                fmt::format("EpochSeries::Since: reference epoch is in {} but the series is in {}",
                            time_to_string.at(ref.scale()), time_to_string.at(scale_)),
                "EpochSeries");
    VecX out(size());
    for (int i = 0; i < size(); i++) {
      double whole = static_cast<double>(sec_[i] - ref.seconds());
      out(i) = Real(whole) + (frac_(i) - ref.fraction());
    }
    return out;
  }

  EpochSeries EpochSeries::To(Time target) const {
    if (target == scale_) return *this;
    std::vector<int64_t> sec = sec_;
    VecX frac(size());
    for (int i = 0; i < size(); i++) {
      Epoch e(sec_[i], frac_(i), scale_);
      frac(i) = frac_(i) + TimeScaleOffset(e, target);
    }
    return EpochSeries(std::move(sec), std::move(frac), target);
  }

}  // namespace lupnt
