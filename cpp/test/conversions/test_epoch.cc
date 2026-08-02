// Epoch: split (int64 seconds + fractional second) representation.
//
// The point of the class is that neither differences nor time-scale
// conversions are limited by the ~0.25 us ULP of an absolute seconds-from-J2000
// double at present-day dates. These tests pin that down.

#include <lupnt/conversions/epoch.h>
#include <lupnt/conversions/time_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // ULP of an absolute epoch [s] at magnitude |t|.
  double EpochUlp(double t) { return std::abs(t) * 2.220446049250313e-16; }
}  // namespace

TEST_CASE("conversions.epoch") {
  // A 2030 epoch: |t| ~ 9.5e8 s, so one ULP is ~0.21 us.
  const Real t_2030 = GregorianToTime("2030-01-01T00:00:00");
  const double ulp = EpochUlp(t_2030.val());

  SECTION("normalisation keeps the fraction in [0,1) without losing precision") {
    Epoch e(1000, Real(2.75), Time::TAI);
    REQUIRE(e.seconds() == 1002);
    REQUIRE_THAT(e.fraction().val(), WithinAbs(0.75, 1.0e-15));

    Epoch neg(1000, Real(-0.25), Time::TAI);
    REQUIRE(neg.seconds() == 999);
    REQUIRE_THAT(neg.fraction().val(), WithinAbs(0.75, 1.0e-15));
  }

  SECTION("differencing two nearby 2030 epochs is exact far below the epoch ULP") {
    // Separate two epochs by 1 picosecond. A bare `Real` cannot represent that
    // separation at all at 2030 (one ULP is ~2e-7 s); Epoch can.
    const double dt = 1.0e-12;
    Epoch a = Epoch::FromSeconds(t_2030, Time::TAI);
    Epoch b = a + Real(dt);

    INFO("absolute-epoch ULP at 2030 = " << ulp << " s");
    REQUIRE_THAT((b - a).val(), WithinRel(dt, 1.0e-9));

    // The same thing done on bare Reals collapses to zero.
    Real a_flat = a.ToSeconds();
    Real b_flat = a_flat + dt;
    REQUIRE((b_flat - a_flat).val() == 0.0);
  }

  SECTION("adding many small steps does not drift") {
    // 1e6 steps of 1e-6 s should land exactly 1 s later.
    Epoch e = Epoch::FromSeconds(t_2030, Time::TAI);
    Epoch start = e;
    for (int i = 0; i < 1000000; i++) e += Real(1.0e-6);
    REQUIRE_THAT((e - start).val(), WithinAbs(1.0, 1.0e-9));
  }

  SECTION("round-tripping through every supported scale returns the same instant") {
    Epoch tai = Epoch::FromSeconds(t_2030, Time::TAI);
    for (Time scale : {Time::UTC, Time::UT1, Time::GPS, Time::TT, Time::TDB, Time::TCG, Time::TCB,
                       Time::TCL, Time::LT}) {
      Epoch there = tai.To(scale);
      Epoch back = there.To(Time::TAI);
      INFO("scale = " << time_to_string.at(scale));
      // Sub-nanosecond round trip, well under the 0.21 us epoch ULP.
      REQUIRE_THAT((back - tai).val(), WithinAbs(0.0, 1.0e-9));
    }
  }

  SECTION("scale conversions agree with the Real-based ConvertTime to within its own ULP") {
    Epoch tai = Epoch::FromSeconds(t_2030, Time::TAI);
    for (Time scale : {Time::UTC, Time::GPS, Time::TT, Time::TDB}) {
      Real reference = ConvertTime(t_2030, Time::TAI, scale);
      Real from_epoch = tai.To(scale).ToSeconds();
      INFO("scale = " << time_to_string.at(scale));
      // Both are absolute epochs here, so they can only agree to a few ULP --
      // that is the floor Epoch exists to avoid, not a defect of either.
      REQUIRE_THAT(from_epoch.val(), WithinAbs(reference.val(), 8.0 * ulp));
    }
  }

  SECTION("TT - TAI is exactly 32.184 s, resolved far below the epoch ULP") {
    Epoch tai = Epoch::FromSeconds(t_2030, Time::TAI);
    Epoch tt = tai.To(Time::TT);
    // Compare in TAI so both sides share a scale.
    Real diff = tt.To(Time::TAI) - tai;
    REQUIRE_THAT(diff.val(), WithinAbs(0.0, 1.0e-9));
    REQUIRE_THAT(TimeScaleOffset(tai, Time::TT).val(), WithinAbs(32.184, 1.0e-12));
    REQUIRE_THAT(TimeScaleOffset(tai, Time::GPS).val(), WithinAbs(-19.0, 1.0e-12));
  }

  SECTION("mixing time scales is rejected rather than silently wrong") {
    Epoch tai = Epoch::FromSeconds(t_2030, Time::TAI);
    Epoch tt = tai.To(Time::TT);
    REQUIRE_THROWS_AS(tt - tai, std::runtime_error);
  }

  SECTION("ordering is exact across a sub-ULP separation") {
    Epoch a = Epoch::FromSeconds(t_2030, Time::TAI);
    Epoch b = a + Real(1.0e-12);
    REQUIRE(a < b);
    REQUIRE(b > a);
    REQUIRE(a != b);
    REQUIRE(a == a);
  }

  SECTION("FromGregorian keeps sub-second precision in the seconds field") {
    Epoch e = Epoch::FromGregorian(2030, 1, 1, 0, 0, Real(30.25), Time::TAI);
    Epoch base = Epoch::FromGregorian(2030, 1, 1, 0, 0, Real(30.0), Time::TAI);
    REQUIRE_THAT((e - base).val(), WithinAbs(0.25, 1.0e-12));
  }
}

TEST_CASE("conversions.epoch_series") {
  const Epoch start = Epoch::FromGregorian("2030-01-01T00:00:00", Time::TAI);

  SECTION("Linspace does not drift over a long grid") {
    // 1e5 samples 0.1 s apart. Building this as t0 + i*dt in a single double
    // at 2030 would quantise every sample to the ~0.21 us epoch grid.
    const int n = 100000;
    EpochSeries s = EpochSeries::Linspace(start, Real(0.1), n);
    REQUIRE(s.size() == n);

    // The last sample must be exactly (n-1)*0.1 s after the first.
    Real elapsed = s[n - 1] - s[0];
    REQUIRE_THAT(elapsed.val(), WithinAbs((n - 1) * 0.1, 1.0e-9));

    // Every sample sits on its exact grid point.
    VecX since = s.Since(start);
    for (int i : {1, 7, 1234, n / 2, n - 1}) {
      INFO("i = " << i);
      REQUIRE_THAT(since(i).val(), WithinAbs(i * 0.1, 1.0e-9));
    }
  }

  SECTION("Since removes the large common epoch exactly") {
    EpochSeries s = EpochSeries::Linspace(start, Real(1.0e-6), 5);
    VecX since = s.Since(start);
    for (int i = 0; i < 5; i++) REQUIRE_THAT(since(i).val(), WithinAbs(i * 1.0e-6, 1.0e-15));
  }

  SECTION("series scale conversion matches the scalar path") {
    EpochSeries s = EpochSeries::Linspace(start, Real(3600.0), 4);
    EpochSeries tdb = s.To(Time::TDB);
    REQUIRE(tdb.scale() == Time::TDB);
    for (int i = 0; i < s.size(); i++) {
      Epoch scalar = s[i].To(Time::TDB);
      INFO("i = " << i);
      REQUIRE_THAT((tdb[i] - scalar).val(), WithinAbs(0.0, 1.0e-12));
    }
  }

  SECTION("round trip through seconds preserves the instants to the epoch ULP") {
    EpochSeries s = EpochSeries::Linspace(start, Real(60.0), 10);
    EpochSeries rt = EpochSeries::FromSeconds(s.ToSeconds(), s.scale());
    for (int i = 0; i < s.size(); i++) {
      REQUIRE_THAT((rt[i] - s[i]).val(), WithinAbs(0.0, 1.0e-6));
    }
  }
}
