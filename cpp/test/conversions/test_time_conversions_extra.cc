#include <lupnt/conversions/time_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <string>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  const double eps_t = 1e-6;
}

// ---------------------------------------------------------------------------
// Fixed, definitional offsets between the atomic/dynamical time scales.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_extra.fixed_offsets") {
  // A handful of epochs [s since J2000]; the offsets are epoch-independent.
  for (double t : {0.0, 1.0e6, -3.5e7}) {
    INFO("t = " << t);
    Real ts(t);

    SECTION("TT - TAI == 32.184 s") {
      REQUIRE_THAT((TaiToTt(ts) - ts).val(), WithinAbs(32.184, 1e-12));
      REQUIRE_THAT((ts - TtToTai(ts)).val(), WithinAbs(32.184, 1e-12));
    }
    SECTION("GPS - TAI == -19 s") {
      REQUIRE_THAT((TaiToGps(ts) - ts).val(), WithinAbs(-19.0, 1e-12));
      REQUIRE_THAT((GpsToTai(ts) - ts).val(), WithinAbs(19.0, 1e-12));
    }
    SECTION("GPS - TT == -51.184 s (TAI-19 then +32.184)") {
      Real t_tt = TaiToTt(ts);
      Real t_gps = TaiToGps(ts);
      REQUIRE_THAT((t_gps - t_tt).val(), WithinAbs(-51.184, 1e-12));
    }
  }
}

// ---------------------------------------------------------------------------
// Round-trip identities across the time-scale conversion pairs.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_extra.round_trips") {
  Real t = 1.2345e7;  // ~143 days past J2000 [s]

  SECTION("TAI <-> TT") { REQUIRE_THAT(TtToTai(TaiToTt(t)).val(), WithinAbs(t.val(), eps_t)); }
  SECTION("TAI <-> GPS") { REQUIRE_THAT(GpsToTai(TaiToGps(t)).val(), WithinAbs(t.val(), eps_t)); }
  SECTION("TT <-> TCG") { REQUIRE_THAT(TcgToTt(TtToTcg(t)).val(), WithinAbs(t.val(), 1e-6)); }
  SECTION("TT <-> TDB") { REQUIRE_THAT(TDBToTt(TtToTdb(t)).val(), WithinAbs(t.val(), 1e-6)); }
  SECTION("TT <-> TCB") {
    Vec3 x0 = Vec3::Zero();
    REQUIRE_THAT(TcbToTt(TtToTcb(t, x0), x0).val(), WithinAbs(t.val(), 1e-6));
  }
  SECTION("TT->TCB->TDB is consistent with the direct TT->TDB conversion") {
    // Going TT -> TCB -> TDB must reproduce the direct TT -> TDB result.
    Real via_tcb = TcbToTdb(TtToTcb(t));
    Real direct = TtToTdb(t);
    REQUIRE_THAT(via_tcb.val(), WithinAbs(direct.val(), 1e-6));
    // TDB stays within a few ms of TT (only periodic + tiny secular terms).
    REQUIRE_THAT(direct.val(), WithinAbs(t.val(), 0.01));
  }

  SECTION("UTC <-> UT1 round trip (needs bundled EOP; skipped if unavailable)") {
    try {
      Real t_utc = 3.0e6;
      Real rt = Ut1ToUtc(UtcToUt1(t_utc));
      REQUIRE_THAT(rt.val(), WithinAbs(t_utc.val(), 1e-3));
    } catch (const std::exception& e) {
      INFO("UTC<->UT1 skipped: EOP data unavailable (" << e.what() << ")");
      SUCCEED();
    }
  }
}

// ---------------------------------------------------------------------------
// MJD / JD / seconds-from-J2000 helpers against documented anchors.
//   J2000 = 2000-01-01 12:00:00 TT = MJD 51544.5 = JD 2451545.0 = 0 s.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_extra.mjd_jd_anchors") {
  SECTION("t = 0 s maps to the J2000 MJD/JD anchors") {
    REQUIRE_THAT(TimeToMjd(Real(0.0)).val(), WithinAbs(51544.5, 1e-9));
    REQUIRE_THAT(TimeToJd(Real(0.0)).val(), WithinAbs(2451545.0, 1e-9));
    REQUIRE_THAT(MjdToTime(Real(51544.5)).val(), WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(JdToTime(Real(2451545.0)).val(), WithinAbs(0.0, 1e-6));
  }

  SECTION("MJD and JD differ by the fixed 2400000.5 offset") {
    Real t = 8.64e5;  // 10 days past J2000
    REQUIRE_THAT((TimeToJd(t) - TimeToMjd(t)).val(), WithinAbs(2400000.5, 1e-9));
  }

  SECTION("MjdToGregorian resolves the J2000 anchor to 2000-01-01 12:00:00") {
    auto [y, mo, d, h, mi, s] = MjdToGregorian(Real(51544.5));
    REQUIRE(y == 2000);
    REQUIRE(mo == 1);
    REQUIRE(d == 1);
    REQUIRE(h == 12);
    REQUIRE(mi == 0);
    REQUIRE_THAT(s.val(), WithinAbs(0.0, 1e-6));
  }

  SECTION("Gregorian string helpers render the J2000 anchor") {
    std::string s_mjd = MjdToGregorianString(Real(51544.5), 3);
    std::string s_time = TimeToGregorianString(Real(0.0), 3);
    REQUIRE(s_mjd.rfind("2000-01-01T12:00:00", 0) == 0);
    REQUIRE(s_time.rfind("2000-01-01T12:00:00", 0) == 0);
  }

  SECTION("GregorianToTime string overload matches the numeric overload") {
    Real t_str = GregorianToTime("2024-03-15T06:30:15.5");
    Real t_num = GregorianToTime(2024, 3, 15, 6, 30, 15.5);
    REQUIRE_THAT(t_str.val(), WithinAbs(t_num.val(), 1e-6));
  }
}

// ---------------------------------------------------------------------------
// Earth-rotation / sidereal-time helpers.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_extra.sidereal") {
  SECTION("EarthRotationAngle advances by ~omega_earth per second, wrapped to [-pi, pi)") {
    Real era0 = EarthRotationAngle(Real(0.0));
    Real era1 = EarthRotationAngle(Real(1.0));  // +1 second UT1
    // EarthRotationAngle wraps to [-pi, pi).
    REQUIRE(era0.val() >= -M_PI);
    REQUIRE(era0.val() < M_PI);
    // Rate: dERA/dt ~ 1.00273781 * (2pi/86400), ~7.292e-5 rad/s.
    Real rate = era1 - era0;
    REQUIRE_THAT(rate.val(), WithinRel(7.2921151467e-5, 1e-3));
  }

  SECTION("GMST and GAST agree to within the equation of the equinoxes (< 1.2e-4 rad)") {
    Real mjd_ut1 = 51544.5 + 10.0;  // 10 days past J2000
    Real gmst = GreenwichMeanSiderealTime(mjd_ut1);
    Real gast = GreenwichApparentSiderealTime(mjd_ut1);
    REQUIRE(gmst.val() >= 0.0);
    REQUIRE(gmst.val() < TWO_PI);
    REQUIRE(gast.val() >= 0.0);
    REQUIRE(gast.val() < TWO_PI);
    // The equation of the equinoxes is at most ~1.1e-4 rad (~1.2 s).
    Real diff = gast - gmst;
    if (diff.val() > M_PI) diff -= TWO_PI;
    if (diff.val() < -M_PI) diff += TWO_PI;
    REQUIRE(std::abs(diff.val()) < 1.2e-4);
  }
}

// ---------------------------------------------------------------------------
// Vectorized overloads must agree elementwise with the scalar functions.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.time_conversions_extra.vectorized") {
  VecX t(3);
  t << 0.0, 1.0e6, -2.5e6;

  SECTION("TaiToTt vectorized matches scalar") {
    VecX out = TaiToTt(t);
    REQUIRE(out.size() == 3);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(out(i).val(), WithinAbs(TaiToTt(Real(t(i))).val(), 1e-9));
  }

  SECTION("MjdToTime vectorized matches scalar") {
    VecX mjd(3);
    mjd << 51544.5, 60000.0, 40000.0;
    VecX out = MjdToTime(mjd);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(out(i).val(), WithinAbs(MjdToTime(Real(mjd(i))).val(), 1e-6));
  }

  SECTION("ConvertTime vectorized matches the scalar ConvertTime") {
    VecX out = ConvertTime(t, Time::TAI, Time::TT);
    REQUIRE(out.size() == 3);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(out(i).val(),
                   WithinAbs(ConvertTime(Real(t(i)), Time::TAI, Time::TT).val(), 1e-9));
  }
}
