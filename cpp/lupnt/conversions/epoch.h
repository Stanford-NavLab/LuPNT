#pragma once
/**
 * @file epoch.h
 * @brief `Epoch` -- a time-system-tagged instant stored as exact integer
 *        seconds plus a fractional second, so differences and time-system
 *        conversions retain sub-nanosecond precision at any date.
 *
 * ## Why this exists
 *
 * LuPNT's default epoch representation is a bare `Real` holding seconds from
 * J2000. At present-day dates |t| ~ 1e9 s, so one double ULP is
 *
 *     |t| * 2^-52  ~  2.5e-7 s  ~  0.25 us   (~75 m at the speed of light)
 *
 * Any small quantity derived from such an epoch -- a time-scale offset, a
 * light time, an integration step -- is therefore snapped to that grid. For a
 * PNT library that is a hard floor on what can be represented, independent of
 * how good the underlying models are.
 *
 * `Epoch` removes the floor by never storing the large and small parts in the
 * same double:
 *
 *     epoch = sec (int64, exact)  +  frac (Real, in [0,1))
 *
 * Differences subtract the integer parts exactly, so `a - b` is a small
 * duration carried at the full precision of a double (~1e-16 s for a
 * one-second separation). Conversions add an inter-scale offset -- always
 * O(100 s) or smaller, hence exact to ~1e-14 s -- to `frac` and renormalise.
 *
 * ## Autodiff contract
 *
 * Derivatives live **entirely on `frac`**; `sec` is a plain integer and is
 * never differentiated. This matches how time enters estimation: the epoch is
 * a fixed parameter, while the small time-like quantities that *are*
 * estimated (clock bias, clock drift) are exactly the ones that stay in
 * `frac`. Renormalisation shifts whole seconds between `sec` and `frac`,
 * which is derivative-preserving because it subtracts an integer constant.
 *
 * ## Interop
 *
 * `ToSeconds()` collapses back to a bare `Real` for the existing API surface.
 * That collapse is lossy by construction (it re-introduces the ULP floor), so
 * prefer `operator-` and `To()` when precision matters.
 */

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"

namespace lupnt {

  /// @brief A time-system-tagged instant: exact integer seconds from J2000
  /// plus a fractional second in [0, 1).
  ///
  /// The time scale is carried with the value so that mixing scales (e.g.
  /// subtracting a TT epoch from a TDB one) is caught rather than silently
  /// producing a wrong answer.
  class Epoch {
  public:
    Epoch() = default;

    /// @brief Construct from integer seconds + fractional second.
    /// The fraction need not be normalised; it is folded into `sec`.
    Epoch(int64_t sec, Real frac, Time scale) : sec_(sec), frac_(frac), scale_(scale) {
      Normalize();
    }

    /// @brief Build from a bare seconds-from-J2000 value.
    ///
    /// Lossy in the sense that `t` has already been rounded to its own ULP;
    /// this cannot recover bits that were never there. Use it at ingestion
    /// (parsing a date, reading a config) rather than mid-computation.
    static Epoch FromSeconds(Real t, Time scale);

    /// @brief Build from a calendar date; the seconds field may be fractional.
    static Epoch FromGregorian(int year, int month, int day, int hour, int min, Real sec,
                               Time scale);

    /// @brief Build from an ISO-like string, e.g. "2025-06-07T12:00:00.000".
    static Epoch FromGregorian(const std::string& date, Time scale);

    int64_t seconds() const { return sec_; }
    Real fraction() const { return frac_; }
    Time scale() const { return scale_; }

    /// @brief Collapse to seconds from J2000.
    ///
    /// **Lossy**: reintroduces the ~0.25 us absolute-epoch ULP at present-day
    /// dates. Provided for interop with the `Real`-based API.
    Real ToSeconds() const { return Real(static_cast<double>(sec_)) + frac_; }

    /// @brief Modified Julian Date in this epoch's scale (lossy, like ToSeconds).
    Real ToMjd() const;

    /// @brief Convert to another time system, retaining sub-nanosecond precision.
    ///
    /// The inter-scale offset is computed as a small quantity (never larger
    /// than ~100 s across the supported scales) and folded into `frac`, so no
    /// large-magnitude cancellation occurs. Supported scales: UT1, UTC, TAI,
    /// GPS, TT, TDB, TCG, TCB, TCL, LT.
    Epoch To(Time target) const;

    /// @brief Format as a calendar string in this epoch's own scale.
    std::string ToGregorianString(int precision = 6) const;

    // --- arithmetic ------------------------------------------------------
    // Adding a duration keeps the scale; subtracting two epochs yields a
    // duration. Durations are plain `Real` seconds -- they are small, so a
    // double represents them exactly.

    Epoch operator+(Real dt) const { return Epoch(sec_, frac_ + dt, scale_); }
    Epoch operator-(Real dt) const { return Epoch(sec_, frac_ - dt, scale_); }
    Epoch& operator+=(Real dt) {
      frac_ += dt;
      Normalize();
      return *this;
    }
    Epoch& operator-=(Real dt) {
      frac_ -= dt;
      Normalize();
      return *this;
    }

    /// @brief Elapsed time `*this - other` [s], exact for nearby epochs.
    ///
    /// The integer parts cancel exactly in int64 arithmetic, so the result
    /// carries the full precision of a double *at the magnitude of the
    /// difference* rather than at the magnitude of the epochs.
    /// Throws if the two epochs are in different time scales.
    Real operator-(const Epoch& other) const;

    // --- comparison ------------------------------------------------------
    // Ordering compares whole seconds first, so it is exact.
    bool operator==(const Epoch& o) const;
    bool operator!=(const Epoch& o) const { return !(*this == o); }
    bool operator<(const Epoch& o) const;
    bool operator>(const Epoch& o) const { return o < *this; }
    bool operator<=(const Epoch& o) const { return !(o < *this); }
    bool operator>=(const Epoch& o) const { return !(*this < o); }

  private:
    /// Fold whole seconds out of `frac_` so that frac_ lies in [0, 1).
    /// Subtracting an integer preserves the autodiff derivative.
    void Normalize();

    int64_t sec_ = 0;
    Real frac_ = 0.0;
    Time scale_ = Time::TAI;
  };

  /// @brief `dt + epoch`, for symmetry with `epoch + dt`.
  inline Epoch operator+(Real dt, const Epoch& e) { return e + dt; }

  /// @brief A sequence of epochs sharing one time scale.
  ///
  /// Much of LuPNT's API is vectorized over epochs (`Propagate(x0, tfs)`,
  /// `GetBodyPosVel(const VecX&, ...)`, every `VEC_DEF_REAL` time function).
  /// `std::vector<Epoch>` would work but loses the contiguous-array shape those
  /// APIs rely on, so `EpochSeries` stores the split representation columnwise:
  /// an exact `int64` second per sample plus a `VecX` of fractions.
  ///
  /// The common case -- a uniform grid -- is built exactly by `Linspace`,
  /// which accumulates in integer seconds rather than by repeated addition of
  /// a rounded step.
  class EpochSeries {
  public:
    EpochSeries() = default;
    EpochSeries(std::vector<int64_t> sec, VecX frac, Time scale)
        : sec_(std::move(sec)), frac_(std::move(frac)), scale_(scale) {
      LUPNT_CHECK(static_cast<int>(sec_.size()) == static_cast<int>(frac_.size()),
                  "EpochSeries: seconds and fraction arrays must be the same length",
                  "EpochSeries");
      Normalize();
    }

    /// @brief `n` epochs spaced `step` seconds apart starting at `start`.
    ///
    /// The i-th sample is `start + i*step` computed in the split
    /// representation, so a long grid does not accumulate rounding the way
    /// `t0 + i*dt` in a single double does.
    static EpochSeries Linspace(const Epoch& start, Real step, int n);

    /// @brief Wrap a bare seconds-from-J2000 vector (lossy, for interop).
    static EpochSeries FromSeconds(const VecX& t, Time scale);

    int size() const { return static_cast<int>(sec_.size()); }
    Time scale() const { return scale_; }
    Epoch operator[](int i) const { return Epoch(sec_[i], frac_(i), scale_); }

    /// @brief Collapse to bare seconds from J2000 (lossy; reintroduces the ULP floor).
    VecX ToSeconds() const;

    /// @brief Elapsed seconds of every sample relative to `ref`, exactly.
    ///
    /// This is the precision-preserving way to feed a vectorized API: the
    /// large common epoch is removed once, and what the callee sees is a set
    /// of small offsets.
    VecX Since(const Epoch& ref) const;

    /// @brief Convert the whole series to another time scale.
    EpochSeries To(Time target) const;

  private:
    void Normalize();

    std::vector<int64_t> sec_;
    VecX frac_;
    Time scale_ = Time::TAI;
  };

  /// @brief Offset `(to - from)` [s] at the instant `e`, as a small quantity.
  ///
  /// This is the primitive `Epoch::To` is built on; exposed because it is also
  /// useful directly. Magnitudes stay below ~100 s for all supported scales,
  /// so the result is exact to ~1e-14 s.
  Real TimeScaleOffset(const Epoch& e, Time to);

}  // namespace lupnt
