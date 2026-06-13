#include <lupnt/conversions/time_conversions.h>
#include <lupnt/data/kernels.h>
#include <lupnt/environment/body.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.time_conversions") {
  SECTION("MJD, JD, and seconds from J2000 round trip") {
    const Real t = 123456.789;

    REQUIRE_THAT(MjdToTime(TimeToMjd(t)).val(), WithinAbs(t.val(), 1.0e-6));
    REQUIRE_THAT(JdToTime(TimeToJd(t)).val(), WithinAbs(t.val(), 2.0e-5));
  }

  SECTION("Gregorian conversion round trips calendar fields") {
    Real mjd = GregorianToMjd(2024, 2, 29, 12, 34, 56.5);
    auto [year, month, day, hour, minute, second] = MjdToGregorian(mjd);

    REQUIRE(year == 2024);
    REQUIRE(month == 2);
    REQUIRE(day == 29);
    REQUIRE(hour == 12);
    REQUIRE(minute == 34);
    REQUIRE_THAT(second.val(), WithinAbs(56.5, 1.0e-6));
  }

  SECTION("TAI, TT, UTC, and GPS conversions are internally consistent") {
    Real t_tai = GregorianToTime(2024, 1, 1, 0, 0, 0.0);

    REQUIRE_THAT(TtToTai(TaiToTt(t_tai)).val(), WithinAbs(t_tai.val(), epsilon));
    REQUIRE_THAT(UtcToTai(TaiToUtc(t_tai)).val(), WithinAbs(t_tai.val(), epsilon));
    REQUIRE_THAT(GpsToTai(TaiToGps(t_tai)).val(), WithinAbs(t_tai.val(), epsilon));
    REQUIRE_THAT(ConvertTime(t_tai, Time::TAI, Time::TAI).val(), WithinAbs(t_tai.val(), epsilon));
  }

  SECTION("coordinate time conversion accepts only coordinate-time scales") {
    Real t_tdb = ConvertTime(GregorianToTime(2024, 1, 1, 0, 0, 0.0), Time::UTC, Time::TDB);

    REQUIRE(IsCoordinateTimeScale(Time::TDB));
    REQUIRE(IsCoordinateTimeScale(Time::TCL));
    REQUIRE_FALSE(IsCoordinateTimeScale(Time::TAI));

    Real t_tcb = ConvertCoordinateTime(t_tdb, Time::TDB, Time::TCB);
    Real t_tdb_roundtrip = ConvertCoordinateTime(t_tcb, Time::TCB, Time::TDB);
    REQUIRE_THAT(t_tdb_roundtrip.val(), WithinAbs(t_tdb.val(), 1.0e-3));

    VecX t_vec(2);
    t_vec << t_tdb, t_tdb + 10.0;
    VecX t_vec_tcb = ConvertCoordinateTime(t_vec, Time::TDB, Time::TCB);
    REQUIRE(t_vec_tcb.size() == 2);
    REQUIRE_THAT(ConvertCoordinateTime(t_vec_tcb(0), Time::TCB, Time::TDB).val(),
                 WithinAbs(t_tdb.val(), 1.0e-3));

    REQUIRE_THROWS_AS(ConvertCoordinateTime(t_tdb, Time::TDB, Time::TAI), std::runtime_error);
  }

  SECTION("coordinate time conversion can route between TDB and TCL") {
    Real t0_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    Real t0_tdb = ConvertCoordinateTime(t0_tcb, Time::TCB, Time::TDB);
    Real t0_tcl = ConvertCoordinateTime(t0_tdb, Time::TDB, Time::TCL);

    REQUIRE_THAT(t0_tcl.val(), WithinAbs(t0_tcb.val(), 1.0e-3));
  }

  SECTION("TCB to TCL uses lunar center by default") {
    Real t0_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    Vec3 x_moon = GetBodyPos(TcbToTdb(t0_tcb), BodyId::SSB, BodyId::MOON, Frame::GCRF);

    REQUIRE_THAT(TcbToTcl(t0_tcb).val(), WithinAbs(TcbToTcl(t0_tcb, x_moon).val(), 1.0e-12));
  }

  SECTION("TCB to TCL position overload uses r_M = x - x_M in BCRS") {
    Real t0_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    Real t0_tdb = TcbToTdb(t0_tcb);
    Vec6 rv_moon = GetBodyPosVel(t0_tdb, BodyId::SSB, BodyId::MOON, Frame::GCRF);
    Vec3 x_moon = rv_moon.head<3>();
    Vec3 v_moon = rv_moon.tail<3>();
    Vec3 x_clock = x_moon + Vec3(1.0e9, 2.0e9, -5.0e8);

    const std::array<BodyId, 9> external_bodies = {
        BodyId::SUN,     BodyId::MERCURY, BodyId::VENUS,  BodyId::EARTH,   BodyId::MARS,
        BodyId::JUPITER, BodyId::SATURN,  BodyId::URANUS, BodyId::NEPTUNE,
    };
    Real w_ext = 0.0;
    for (BodyId body : external_bodies) {
      Vec3 x_body = GetBodyPos(t0_tdb, BodyId::SSB, body, Frame::GCRF);
      w_ext += GetBodyGM(body) / (x_moon - x_body).norm();
    }

    Real v_dot_r = v_moon.dot(x_clock - x_moon);
    Real c2 = C * C;
    Real expected
        = -v_dot_r / c2 - (0.5 * v_moon.squaredNorm() + 3.0 * w_ext) * v_dot_r / (c2 * c2);

    Real actual = TcbToTcl(t0_tcb, x_clock) - TcbToTcl(t0_tcb);
    REQUIRE_THAT(actual.val(), WithinAbs(expected.val(), 1.0e-7));
  }

  SECTION("TCB and TCL round trip with and without BCRS position") {
    Real t_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB) + 2.0 * SECS_DAY;
    Vec6 rv_moon = GetBodyPosVel(TcbToTdb(t_tcb), BodyId::SSB, BodyId::MOON, Frame::GCRF);
    Vec3 x_clock = rv_moon.head<3>() + Vec3(100.0, -250.0, 75.0);

    REQUIRE_THAT(TclToTcb(TcbToTcl(t_tcb)).val(), WithinAbs(t_tcb.val(), 1.0e-9));
    REQUIRE_THAT(TclToTcb(TcbToTcl(t_tcb, x_clock), x_clock).val(), WithinAbs(t_tcb.val(), 1.0e-9));

    VecX t_vec(2);
    t_vec << t_tcb, t_tcb + 10.0;
    MatX3 x_vec(2, 3);
    x_vec.row(0) = x_clock.transpose();
    x_vec.row(1) = (x_clock + Vec3(1.0, 2.0, 3.0)).transpose();
    VecX t_tcl = ConvertCoordinateTime(t_vec, Time::TCB, Time::TCL, x_vec);
    VecX t_rt = ConvertCoordinateTime(t_tcl, Time::TCL, Time::TCB, x_vec);
    REQUIRE_THAT(t_rt(0).val(), WithinAbs(t_vec(0).val(), 1.0e-9));
    REQUIRE_THAT(t_rt(1).val(), WithinAbs(t_vec(1).val(), 1.0e-9));
  }
}

// ============================================================================
// TL − TT tests (Turyshev 2026, Eq. 57)
//
// All tests use t_tdb near T_0 so the numerical integration only spans a few
// days, keeping test run-times fast.  T_0 ≈ -725,803,194 s from J2000.
// ============================================================================

TEST_CASE("conversions.lt_minus_tt") {
  // Reference epoch: coordinate time origin T_0 (1977-01-01 00:00:32.184 TAI)
  const Real t0_tdb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);  // ≈ -7.258e8 s
  const Real t0_eq57 = t0_tdb + TDB_0;
  auto earth_moon_endpoint = [](Real t_tdb) {
    Vec6 rv_em = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec6 rv_earth = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
    return rv_earth.tail<3>().dot(rv_em.head<3>());
  };

  SECTION("TL-TT is exactly zero at the reference epoch T_0") {
    // At T_0, the integration interval is empty → TL-TT = secular × 0 = 0
    Real lt_minus_tt = TdbToLtMinusTt(t0_tdb);
    REQUIRE_THAT(lt_minus_tt.val(), WithinAbs(0.0, 1.0e-12));
  }

  SECTION("TL-TT secular rate: net rate is (L_G-L_L)/(1-L_B) - L_EM over 1 day") {
    // The net secular rate is (L_G-L_L)/(1-L_B) minus the mean c^{-2} integrand
    // (≈ L_EM = 1.709e-11 s/s), plus the Eq. 56 endpoint contribution.
    const double dt = 1.0 * SECS_DAY;
    Real t1 = t0_tdb + dt;
    Real lt_minus_tt_1 = TdbToLtMinusTt(t1);

    double net_rate = (L_G - L_L) / (1.0 - L_B) - L_EM;
    double endpoint = ((earth_moon_endpoint(t0_eq57) - earth_moon_endpoint(t1)) / (C * C)).val();
    double expected_net = net_rate * (t1 - t0_eq57).val() + endpoint;

    REQUIRE_THAT(lt_minus_tt_1.val(), WithinAbs(expected_net, 2.0e-6));
  }

  SECTION("TL-TT net secular rate holds over 30 days") {
    // The endpoint term is periodic, so include it explicitly when checking the
    // underlying secular rate over a longer arc.
    const double dt = 30.0 * SECS_DAY;
    Real t1 = t0_tdb + dt;
    Real lt_minus_tt = TdbToLtMinusTt(t1);

    double net_rate = (L_G - L_L) / (1.0 - L_B) - L_EM;
    double endpoint = ((earth_moon_endpoint(t0_eq57) - earth_moon_endpoint(t1)) / (C * C)).val();
    double expected_net = net_rate * (t1 - t0_eq57).val() + endpoint;

    REQUIRE_THAT(lt_minus_tt.val(), WithinAbs(expected_net, 2.0e-6));
  }

  SECTION("ConvertTime TDB→LT matches TdbToLt exactly") {
    // ConvertTime must route TDB→LT through TdbToLt; results should be bitwise equal.
    const double dt = 30.0 * SECS_DAY;
    Real t_tdb = t0_tdb + dt;

    Real tl_direct = TdbToLt(t_tdb);
    Real tl_route = ConvertTime(t_tdb, Time::TDB, Time::LT);

    REQUIRE_THAT(tl_route.val(), WithinAbs(tl_direct.val(), 1.0e-12));
  }

  SECTION("TDB to TCL route uses TCB to TCL conversion") {
    const double dt = 30.0 * SECS_DAY;
    Real t_tdb = t0_tdb + dt;

    Real t_tcb = ConvertTime(t_tdb, Time::TDB, Time::TCB);
    Real t_tcl_direct = TcbToTcl(t_tcb);
    Real t_tcl_route = ConvertTime(t_tdb, Time::TDB, Time::TCL);

    REQUIRE_THAT(t_tcl_route.val(), WithinAbs(t_tcl_direct.val(), 1.0e-12));
  }

  SECTION("TT and TL agree through Eq.57 and TCB/TCL paths") {
    for (double days : {1.0, 7.0, 30.0}) {
      Real t_tt = t0_tdb + days * SECS_DAY;
      Real t_tdb_guess = TtToTdb(t_tt);
      Vec3 x_earth = GetBodyPos(t_tdb_guess, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Real t_tdb = TtToTdb(t_tt, x_earth);
      x_earth = GetBodyPos(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Vec3 x_moon = GetBodyPos(t_tdb_guess, BodyId::SSB, BodyId::MOON, Frame::GCRF);
      x_moon = GetBodyPos(t_tdb, BodyId::SSB, BodyId::MOON, Frame::GCRF);

      Real t_lt_eq57 = t_tt + TdbToLtMinusTt(t_tdb);
      Real t_tcb = TtToTcb(t_tt, x_earth);
      Real t_lt_tcb = TclToLt(TcbToTcl(t_tcb, x_moon));

      INFO("days = " << days);
      REQUIRE_THAT(t_lt_tcb.val(), WithinAbs(t_lt_eq57.val(), 1.0e-9));

      Real t_tdb_from_lt = t_tdb;
      for (int iter = 0; iter < 10; ++iter) {
        Vec3 x_earth_iter = GetBodyPos(t_tdb_from_lt, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
        Real t_lt_est = TDBToTt(t_tdb_from_lt, x_earth_iter) + TdbToLtMinusTt(t_tdb_from_lt);
        Real delta = t_lt_eq57 - t_lt_est;
        t_tdb_from_lt += delta;
        if (std::abs(delta.val()) < 1.0e-12) break;
      }
      Vec3 x_earth_rt = GetBodyPos(t_tdb_from_lt, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Real t_tt_eq57 = TDBToTt(t_tdb_from_lt, x_earth_rt);
      Real t_tt_tcb = TcbToTt(TclToTcb(LtToTcl(t_lt_eq57), x_moon), x_earth);

      REQUIRE_THAT(t_tt_eq57.val(), WithinAbs(t_tt.val(), 1.0e-9));
      REQUIRE_THAT(t_tt_tcb.val(), WithinAbs(t_tt.val(), 1.0e-9));
    }
  }

  SECTION("TDB → LT → TDB round-trip to within 1 ns") {
    const double dt = 30.0 * SECS_DAY;
    Real t_tdb = t0_tdb + dt;
    Real t_lt = TdbToLt(t_tdb);
    Real t_tdb_rt = LtToTdb(t_lt);
    REQUIRE_THAT(t_tdb_rt.val(), WithinAbs(t_tdb.val(), 1.0e-9));
  }

  SECTION("Chebyshev fit matches direct integration to within 1 ns") {
    // Fit a 3-day window starting at T_0 with 1-day segments, 13 coefficients
    const double dt_win = 3.0 * SECS_DAY;
    Real t_start = t0_tdb;
    Real t_end = t0_tdb + dt_win;

    InitLtMinusTtFit(t_start, t_end, SECS_DAY, 13);

    // Sample several epochs inside the window
    for (int i = 1; i <= 5; i++) {
      Real t_tdb = t0_tdb + i * 0.5 * SECS_DAY;
      REQUIRE(HasFittedLtMinusTt(t_tdb));

      Real direct = TdbToLtMinusTt(t_tdb);  // uses the fit (fast path)
      // Direct integration reference: temporarily clear the fit
      ClearLtMinusTtFit();
      Real reference = TdbToLtMinusTt(t_tdb);
      // Restore fit for remaining iterations
      InitLtMinusTtFit(t_start, t_end, SECS_DAY, 13);

      REQUIRE_THAT(direct.val(), WithinAbs(reference.val(), 1.0e-9));
    }

    ClearLtMinusTtFit();
    REQUIRE_FALSE(HasFittedLtMinusTt(t0_tdb + 1.0 * SECS_DAY));
  }
}
