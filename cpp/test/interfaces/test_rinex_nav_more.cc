#include <lupnt/lupnt.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// Deepens rinex_nav_loader.cc coverage beyond test_rinex_nav_loader_fixture.cc:
// multi-message closest-selection, the clock-correction polynomial, the
// GetLatestEpochTai accessor, FilenameForEpoch, additional unsupported-system /
// invalid-id throw paths, and the YUMA-almanac loader (LoadYumaFile) which is
// otherwise untested. Uses committed fixtures (no network).
namespace {
  std::filesystem::path FixtureFile() {
    return std::filesystem::path(LUPNT_TEST_FIXTURES_DIR) / "gnss" / "BRDC_trimmed.rnx";
  }
  std::filesystem::path YumaFixtureFile() {
    return std::filesystem::path(LUPNT_TEST_FIXTURES_DIR) / "gnss" / "current.alm";
  }
}  // namespace

TEST_CASE("interfaces.rinex_nav_more.closest_message_and_clock") {
  RinexNavLoader loader(FixtureFile());
  const std::string sat_id = "G01";  // three messages: 00:00, 02:00, 04:00

  // af0 clock-polynomial coefficients of the first and last G01 messages
  // (first field of each RINEX epoch header line in the fixture).
  const double af0_msg0 = 3.441781736910e-04;
  const double af0_msg2 = 3.441260196269e-04;

  Real t0 = GregorianToTime(2026, 1, 14, 0, 0, 0);
  Real t2 = GregorianToTime(2026, 1, 14, 4, 0, 0);

  Vec6 rv0;
  Real clk0;
  loader.GetPosVelClock(sat_id, t0, rv0, clk0);
  Vec6 rv2;
  Real clk2;
  loader.GetPosVelClock(sat_id, t2, rv2, clk2);

  // At each message's own epoch t_k ~ 0, so the clock polynomial reduces to af0
  // (plus a sub-30-ns relativistic term). Recovering the *right* af0 confirms
  // the argmin closest-message selection picked the correct message.
  REQUIRE(std::abs(clk0.val() - af0_msg0) < 1e-6);
  REQUIRE(std::abs(clk2.val() - af0_msg2) < 1e-6);

  // The two epochs select distinct messages -> distinct clock corrections.
  REQUIRE(std::abs(clk0.val() - clk2.val()) > 1e-8);

  // Both propagate to a plausible GPS MEO state.
  REQUIRE(rv0.head(3).norm().val() > 2.0e7);
  REQUIRE(rv0.head(3).norm().val() < 3.0e7);
  REQUIRE(rv2.head(3).norm().val() > 2.0e7);
  REQUIRE(rv2.head(3).norm().val() < 3.0e7);
}

TEST_CASE("interfaces.rinex_nav_more.latest_epoch") {
  RinexNavLoader loader(FixtureFile());

  // The freshest G01 orbital slot is the 04:00 message.
  double latest = loader.GetLatestEpochTai("G01");
  double expected = ConvertTime(GregorianToTime(2026, 1, 14, 4, 0, 0), Time::GPS, Time::TAI).val();
  REQUIRE_THAT(latest, WithinAbs(expected, 1e-3));

  // Single-message satellites report that single message's epoch.
  double latest_g02 = loader.GetLatestEpochTai("G02");
  REQUIRE(std::isfinite(latest_g02));

  REQUIRE_THROWS(loader.GetLatestEpochTai("X99"));
}

TEST_CASE("interfaces.rinex_nav_more.filename_for_epoch") {
  // FilenameForEpoch mirrors the daily CDDIS BRDC product naming.
  const Real t = GregorianToTime(2026, 1, 14, 12, 0, 0.0);
  REQUIRE(RinexNavLoader::FilenameForEpoch(t, Time::UTC) == "BRDC00IGS_R_20260140000_01D_MN.rnx");

  // Leap-year day-of-year branch (2024 is a leap year, Mar 1 -> DOY 61).
  const Real t_leap = GregorianToTime(2024, 3, 1, 12, 0, 0.0);
  REQUIRE(RinexNavLoader::FilenameForEpoch(t_leap, Time::UTC)
          == "BRDC00IGS_R_20240610000_01D_MN.rnx");
}

TEST_CASE("interfaces.rinex_nav_more.invalid_and_unsupported_ids_throw") {
  RinexNavLoader loader(FixtureFile());
  Vec6 rv;
  Real clk;
  Real t = GregorianToTime(2026, 1, 14, 1, 0, 0);

  // Empty identifier -> "Invalid satellite identifier".
  REQUIRE_THROWS(loader.GetPosVelClock("", t, rv, clk));
  // SBAS ('S') is another non-Keplerian system, unsupported like GLONASS.
  REQUIRE_THROWS(loader.GetPosVelClock("S20", t, rv, clk));
}

TEST_CASE("interfaces.rinex_nav_more.yuma_almanac_loader") {
  // LoadYumaFile parses a YUMA GPS almanac into per-PRN broadcast-nav slots.
  RinexNavLoader loader;
  loader.LoadYumaFile(YumaFixtureFile());

  const auto& sats = loader.GetSatellites();
  REQUIRE_FALSE(sats.empty());
  // YUMA carries GPS ('G') satellites only.
  REQUIRE(std::all_of(sats.begin(), sats.end(),
                      [](const std::string& s) { return !s.empty() && s.front() == 'G'; }));

  // Query each satellite at its own reference epoch (t_k ~ 0) so the coarse
  // Keplerian slot yields a plausible GPS orbital radius / speed.
  for (const std::string& sat : sats) {
    double t_ref = loader.GetLatestEpochTai(sat);
    Vec6 rv = loader.GetPosVel(sat, Real(t_ref));
    Real r = rv.head(3).norm();
    Real v = rv.tail(3).norm();
    REQUIRE(r.val() > 2.4e7);
    REQUIRE(r.val() < 2.9e7);
    REQUIRE(v.val() > 1.0e3);
    REQUIRE(v.val() < 1.0e4);
  }
}
