#include <lupnt/conversions/time_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <string>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ---------------------------------------------------------------------------
// ConvertTime dispatcher: identity and multi-hop routing through the time-scale
// graph. Only the atomic/dynamical scales (no UTC/UT1, which need EOP; no
// TCL/LT, which need ephemerides) are exercised here.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_more.convert_time_identity") {
  Real t = 4.56e6;
  for (Time scale : {Time::TAI, Time::TT, Time::GPS, Time::TDB, Time::TCG, Time::TCB}) {
    INFO("scale = " << static_cast<int>(scale));
    REQUIRE_THAT(ConvertTime(t, scale, scale).val(), WithinAbs(t.val(), 1e-12));
  }
}

TEST_CASE("conversions.time_conversions_more.convert_time_matches_direct_chains") {
  Real t = 2.5e7;

  // ConvertTime now delegates to Epoch, which composes small offsets and only
  // materialises an absolute epoch at the end; the hand-rolled chains below form
  // an absolute at every step. Both are correct, but they round differently, so
  // they can only agree to the ABSOLUTE-EPOCH ULP -- 5.5 ns at t = 2.5e7 s. The
  // old 1e-9 tolerance was below the representable resolution and passed only
  // because the two routes happened to round identically.
  const double kUlp = std::abs(t.val()) * 2.220446049250313e-16;
  const double kTol = std::max(1.0e-9, 4.0 * kUlp);

  SECTION("TT -> TDB equals the direct TtToTdb") {
    REQUIRE_THAT(ConvertTime(t, Time::TT, Time::TDB).val(), WithinAbs(TtToTdb(t).val(), kTol));
  }
  SECTION("TAI -> TDB equals TtToTdb(TaiToTt)") {
    REQUIRE_THAT(ConvertTime(t, Time::TAI, Time::TDB).val(),
                 WithinAbs(TtToTdb(TaiToTt(t)).val(), kTol));
  }
  SECTION("GPS -> TT equals TaiToTt(GpsToTai)") {
    REQUIRE_THAT(ConvertTime(t, Time::GPS, Time::TT).val(),
                 WithinAbs(TaiToTt(GpsToTai(t)).val(), kTol));
  }
  SECTION("TCG -> TAI equals TtToTai(TcgToTt)") {
    REQUIRE_THAT(ConvertTime(t, Time::TCG, Time::TAI).val(),
                 WithinAbs(TtToTai(TcgToTt(t)).val(), kTol));
  }
  SECTION("TCB -> GPS equals routing through TDB/TT/TAI") {
    Real via = TaiToGps(TtToTai(TDBToTt(TcbToTdb(t))));
    REQUIRE_THAT(ConvertTime(t, Time::TCB, Time::GPS).val(), WithinAbs(via.val(), 1e-6));
  }
}

TEST_CASE("conversions.time_conversions_more.convert_time_round_trips") {
  Real t = -8.0e6;
  // Each pairing routes through a different set of switch arms. Coordinate-time
  // round trips accumulate small numerical error through the (iterative) TCB/TDB
  // inversions, so a millisecond-level tolerance is used.
  REQUIRE_THAT(ConvertTime(ConvertTime(t, Time::TDB, Time::TCG), Time::TCG, Time::TDB).val(),
               WithinAbs(t.val(), 1e-3));
  REQUIRE_THAT(ConvertTime(ConvertTime(t, Time::GPS, Time::TCB), Time::TCB, Time::GPS).val(),
               WithinAbs(t.val(), 1e-3));
  REQUIRE_THAT(ConvertTime(ConvertTime(t, Time::TCG, Time::TCB), Time::TCB, Time::TCG).val(),
               WithinAbs(t.val(), 1e-3));
}

// ---------------------------------------------------------------------------
// ConvertCoordinateTime: the position-independent 3-arg overload, the
// position-aware 4-arg overload, and their vectorized forms. Also the guard
// that rejects non-coordinate time scales.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_more.convert_coordinate_time") {
  Real t = 1.0e7;
  Vec3 x0 = Vec3::Zero();

  SECTION("3-arg overload matches ConvertTime for coordinate scales") {
    REQUIRE_THAT(ConvertCoordinateTime(t, Time::TT, Time::TDB).val(),
                 WithinAbs(ConvertTime(t, Time::TT, Time::TDB).val(), 1e-12));
  }

  SECTION("4-arg TT<->TDB round trip at zero position") {
    Real tdb = ConvertCoordinateTime(t, Time::TT, Time::TDB, x0);
    Real tt = ConvertCoordinateTime(tdb, Time::TDB, Time::TT, x0);
    REQUIRE_THAT(tt.val(), WithinAbs(t.val(), 1e-6));
  }

  SECTION("4-arg identity returns input") {
    REQUIRE_THAT(ConvertCoordinateTime(t, Time::TT, Time::TT, x0).val(), WithinAbs(t.val(), 1e-12));
  }

  SECTION("non-coordinate scale is rejected") {
    REQUIRE_THROWS(ConvertCoordinateTime(t, Time::TAI, Time::TT));
    REQUIRE_THROWS(ConvertCoordinateTime(t, Time::TT, Time::GPS));
  }

  SECTION("vectorized overloads agree elementwise with the corresponding scalar") {
    VecX ts(3);
    ts << 0.0, 1.0e6, -2.0e6;

    // 3-arg vectorized matches the 3-arg scalar.
    VecX out = ConvertCoordinateTime(ts, Time::TT, Time::TDB);
    REQUIRE(out.size() == 3);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(out(i).val(),
                   WithinAbs(ConvertCoordinateTime(Real(ts(i)), Time::TT, Time::TDB).val(), 1e-9));

    // 4-arg vectorized (zero positions) matches the 4-arg scalar per row.
    MatX3 pos = MatX3::Zero(3, 3);
    VecX out_pos = ConvertCoordinateTime(ts, Time::TT, Time::TDB, pos);
    REQUIRE(out_pos.size() == 3);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(
          out_pos(i).val(),
          WithinAbs(ConvertCoordinateTime(Real(ts(i)), Time::TT, Time::TDB, x0).val(), 1e-9));
  }
}

// ---------------------------------------------------------------------------
// Gregorian <-> MJD calendar helpers.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_more.gregorian_to_mjd") {
  SECTION("J2000 civil epoch maps to MJD 51544.5") {
    REQUIRE_THAT(GregorianToMjd(2000, 1, 1, 12, 0, Real(0.0)).val(), WithinAbs(51544.5, 1e-9));
  }

  SECTION("midnight of an MJD-anchor date has an integer MJD") {
    // 1858-11-17T00:00:00 is MJD 0 by definition.
    REQUIRE_THAT(GregorianToMjd(1858, 11, 17, 0, 0, Real(0.0)).val(), WithinAbs(0.0, 1e-9));
  }

  SECTION("GregorianToMjd and MjdToGregorian invert each other at a generic date") {
    Real mjd = GregorianToMjd(2024, 3, 15, 6, 30, Real(15.0));
    auto [y, mo, d, h, mi, s] = MjdToGregorian(mjd);
    REQUIRE(y == 2024);
    REQUIRE(mo == 3);
    REQUIRE(d == 15);
    REQUIRE(h == 6);
    REQUIRE(mi == 30);
    REQUIRE_THAT(s.val(), WithinAbs(15.0, 1e-3));
  }

  SECTION("half-day fraction advances the MJD by 0.5") {
    Real noon = GregorianToMjd(2020, 6, 1, 12, 0, Real(0.0));
    Real midnight = GregorianToMjd(2020, 6, 1, 0, 0, Real(0.0));
    REQUIRE_THAT((noon - midnight).val(), WithinAbs(0.5, 1e-9));
  }
}

// ---------------------------------------------------------------------------
// Gregorian string rendering: the precision argument controls the number of
// fractional-second digits.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_more.gregorian_string_precision") {
  // A time with a fractional-second component (t is TT seconds past J2000).
  Real t = 0.25;  // 0.25 s past the J2000 civil-noon anchor
  std::string p0 = TimeToGregorianString(t, 0);
  std::string p3 = TimeToGregorianString(t, 3);
  std::string p6 = TimeToGregorianString(t, 6);

  // All render the same date/time prefix.
  REQUIRE(p0.rfind("2000-01-01T12:00:00", 0) == 0);
  REQUIRE(p3.rfind("2000-01-01T12:00:00", 0) == 0);

  // Higher precision produces a longer (more decimal digits) string.
  REQUIRE(p6.size() > p3.size());
  REQUIRE(p3.size() > p0.size());
  // Only the >0-precision renderings carry a fractional-second separator.
  REQUIRE_THAT(p3, ContainsSubstring("."));
  REQUIRE(p0.find('.') == std::string::npos);
}
