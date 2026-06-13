#include <lupnt/conversions/time_conversions.h>
#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation against Orekit
//
// The reference values below were generated once via Orekit and checked
// into `cpp/test/orekit/data/orekit_reference.json` by
// `cpp/test/orekit/gen_orekit_reference.py`. This test does **not**
// require Orekit/Java to run -- see `cpp/test/orekit/README.md` for how
// the fixture was produced and how to regenerate it.
// ============================================================================

TEST_CASE("conversions.time_conversions_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  for (const auto& tc : data["time_scales"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    Real t_utc = GregorianToTime(epoch);

    DYNAMIC_SECTION("epoch " << epoch) {
      Real t_tai = ConvertTime(t_utc, Time::UTC, Time::TAI);
      Real t_tt = ConvertTime(t_utc, Time::UTC, Time::TT);
      Real t_tdb = ConvertTime(t_utc, Time::UTC, Time::TDB);
      Real t_gps = ConvertTime(t_utc, Time::UTC, Time::GPS);
      Real t_ut1 = ConvertTime(t_utc, Time::UTC, Time::UT1);

      // TAI-UTC (leap seconds), TT-TAI (exact 32.184 s), and GPS-TAI (exact
      // -19 s) are constant/integer offsets defined by the time-scale
      // standards themselves, so algorithmically LuPNT and Orekit should
      // agree exactly. In practice, LuPNT represents time as "seconds since
      // J2000" in a double, so e.g. `t_tt - t_tai` is computed as
      // `(t_tai + 32.184) - t_tai`. For epochs far from J2000 (the
      // 2024-03-15 case is ~7.6e8 s, i.e. ~24 years, from J2000), the ULP of
      // `t_tai` is ~1.2e-7 s, so adding 32.184 and subtracting back can be
      // off from the exact constant by up to ~half a ULP (~6e-8 s). This is
      // a floating-point representation artifact, not an algorithmic
      // discrepancy -- 1e-6 s (1 microsecond) comfortably absorbs it for
      // epochs within centuries of J2000 while still catching any real
      // (e.g. wrong leap-second count, wrong sign) error.
      RequireNear(t_tai - t_utc, Real(tc["tai_minus_utc"].get<double>()), 1e-6);
      RequireNear(t_tt - t_tai, Real(tc["tt_minus_tai"].get<double>()), 1e-6);
      RequireNear(t_gps - t_tai, Real(tc["gps_minus_tai"].get<double>()), 1e-6);

      // TDB-TAI has a small (~ms) periodic term that depends on the
      // analytical/ephemeris model used for the relativistic correction;
      // LuPNT and Orekit agree to ~30 ns in practice, so 1 us is generous.
      RequireNear(t_tdb - t_tai, Real(tc["tdb_minus_tai"].get<double>()), 1e-6);

      // UT1-UTC depends on the Earth Orientation Parameters (EOP) table in
      // use. LuPNT's bundled EOP table and Orekit's bundled `orekit-data`
      // table are both valid IERS products but correspond to different
      // bulletin "vintages", so they can disagree by tens of milliseconds
      // for the same epoch -- this is a data-vintage difference, not an
      // algorithmic bug. UT1-UTC is bounded to +/-0.9 s by definition (IERS
      // leap-second policy), so a 0.1 s tolerance is loose enough to absorb
      // table differences while still catching gross sign/scale/unit
      // errors.
      //
      // Epochs flagged `near_leap_second` (within ~1 day of a leap-second
      // insertion) are excluded: LuPNT's EOP table stores daily UT1-UTC
      // values and linearly interpolates *across* the 1 s leap jump, which
      // produces up to ~0.94 s of spurious UT1-UTC one hour before an
      // insertion (Orekit interpolates the continuous UT1-TAI instead and
      // is unaffected). This is a known LuPNT limitation, documented in
      // docs/pages/cross_validation.rst.
      if (!tc["near_leap_second"].get<bool>()) {
        RequireNear(t_ut1 - t_utc, Real(tc["ut1_minus_utc"].get<double>()), 0.1);
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Sidereal / Earth-rotation angles
// ----------------------------------------------------------------------------
TEST_CASE("conversions.sidereal_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  for (const auto& tc : data["sidereal"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    // GMST/ERA are functions of UT1, so epochs flagged `near_leap_second`
    // are skipped for the same reason as the UT1-UTC check above (LuPNT's
    // daily UT1-UTC interpolation across the leap jump).
    if (tc["near_leap_second"].get<bool>()) continue;
    Real t_utc = GregorianToTime(epoch);
    Real t_ut1 = ConvertTime(t_utc, Time::UTC, Time::UT1);

    DYNAMIC_SECTION("epoch " << epoch) {
      // GMST: LuPNT implements the classic IAU-82 polynomial, which is the
      // same expression Orekit uses for the IERS 1996 conventions, so the
      // *formulas* agree exactly. The remaining difference comes from the
      // UT1 input itself: LuPNT's and Orekit's EOP tables are different
      // (but both valid) IERS bulletin vintages, and a UT1-UTC difference
      // `dUT1` shifts the angle by `omega_earth * dUT1`
      // (~7.3e-5 rad/s * tens of ms ~ a few microradians). The same applies
      // to the IAU-2000 Earth Rotation Angle, which is a linear function of
      // UT1 with identical coefficients on both sides. 1e-5 rad covers a
      // ~0.14 s UT1 difference (consistent with the 0.1 s UT1-UTC tolerance
      // above) and still catches wrong-formula/wrong-time-scale errors,
      // which would be off by orders of magnitude more.
      Real mjd_ut1 = TimeToMjd(t_ut1);
      Real gmst = GreenwichMeanSiderealTime(mjd_ut1);
      Real era = EarthRotationAngle(t_ut1);

      Real gmst_ref = tc["gmst"].get<double>();
      Real era_ref = tc["era"].get<double>();
      // Compare modulo 2*pi (the fixture stores [0, 2*pi), LuPNT may wrap
      // to (-pi, pi]).
      RequireNear(WrapToPi(gmst - gmst_ref), Real(0.0), 1e-5);
      RequireNear(WrapToPi(era - era_ref), Real(0.0), 1e-5);
    }
  }
}
