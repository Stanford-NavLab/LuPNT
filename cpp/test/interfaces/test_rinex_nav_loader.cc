#include <lupnt/lupnt.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // See `test_sp3_loader.cc` for why this differs from `GetOutputDir("gnss_files")`.
  std::filesystem::path GnssFilesDir() {
    return GetDataPath().parent_path().parent_path() / "output" / "gnss_files";
  }
}  // namespace

const double epsilon = 1e-6;

TEST_CASE("interfaces.rinex_nav_loader") {
  const std::filesystem::path brdc_dir = GnssFilesDir() / "brdc";
  const std::filesystem::path nav_file = brdc_dir / "BRDC00IGS_R_20260140000_01D_MN.rnx";

  const std::string sat_id = "G01";  // GPS PRN 01: navigation messages start 2026-01-14 00:00

  Real t_tai = GregorianToTime(2026, 1, 14, 1, 0, 0);

  SECTION("Default-constructed loader has no satellites loaded") {
    RinexNavLoader loader;
    REQUIRE(loader.GetSatellites().empty());
    REQUIRE_FALSE(loader.HasSatellite(sat_id));
  }

  SECTION("Single-file constructor parses broadcast navigation messages") {
    RinexNavLoader loader(nav_file);

    REQUIRE_FALSE(loader.GetSatellites().empty());
    REQUIRE(loader.HasSatellite(sat_id));
    REQUIRE_FALSE(loader.HasSatellite("X99"));

    // GLONASS ('R', tabulated state-vector ephemeris) is intentionally excluded
    const auto& sats = loader.GetSatellites();
    REQUIRE(std::none_of(sats.begin(), sats.end(),
                         [](const std::string& s) { return !s.empty() && s.front() == 'R'; }));

    Vec6 rv_ecef;
    Real clock_corr_s;
    loader.GetPosVelClock(sat_id, t_tai, rv_ecef, clock_corr_s);

    // GPS broadcast-ephemeris position/velocity should be near the nominal orbit
    Real r = rv_ecef.head(3).norm();
    REQUIRE(r.val() > 2.0e7);
    REQUIRE(r.val() < 3.0e7);

    Real v = rv_ecef.tail(3).norm();
    REQUIRE(v.val() > 1.0e3);
    REQUIRE(v.val() < 1.0e4);

    // Broadcast clock (polynomial + relativistic) corrections are sub-millisecond
    REQUIRE(std::abs(clock_corr_s.val()) < 1e-3);

    // GetPosVel returns the same position/velocity as GetPosVelClock
    Vec6 rv_only = loader.GetPosVel(sat_id, t_tai);
    RequireNear(rv_only, rv_ecef, epsilon);
  }

  SECTION("LoadFile / multi-file constructor accumulate navigation messages") {
    RinexNavLoader loader_a(nav_file);
    REQUIRE(loader_a.HasSatellite(sat_id));

    // Re-parsing the same file concatenates messages but should not change the
    // broadcast position/velocity at a fixed epoch (closest-message selection).
    loader_a.LoadFile(nav_file);
    REQUIRE(loader_a.HasSatellite(sat_id));
    Vec6 rv_a = loader_a.GetPosVel(sat_id, t_tai);

    RinexNavLoader loader_b(std::vector<std::filesystem::path>{nav_file});
    REQUIRE(loader_b.HasSatellite(sat_id));
    Vec6 rv_b = loader_b.GetPosVel(sat_id, t_tai);

    RequireNear(rv_a, rv_b, epsilon);
  }

  SECTION("Queries for unknown or unsupported (GLONASS) satellites throw") {
    RinexNavLoader loader(nav_file);
    Vec6 rv_ecef;
    Real clock_corr_s;

    REQUIRE_THROWS(loader.GetPosVelClock("X99", t_tai, rv_ecef, clock_corr_s));
    // GLONASS ('R'): tabulated state vectors, unsupported by Keplerian propagation
    REQUIRE_THROWS(loader.GetPosVelClock("R01", t_tai, rv_ecef, clock_corr_s));
  }
}
