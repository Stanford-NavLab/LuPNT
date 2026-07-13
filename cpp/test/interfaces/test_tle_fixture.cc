#include <lupnt/interfaces/tle.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>

#include "../utils.cc"
#include "lupnt/agents/agent.h"
#include "lupnt/core/file.h"

using namespace lupnt;
using namespace Catch::Matchers;

// Exercises the interfaces/tle.cc parser (TLE::FromLines / GetPrn), which is
// otherwise only reached from the CI-excluded states.tle test. Inline TLE
// strings keep this data-independent; a bundled-file section (FromFile /
// LoadTleFile) is guarded so it self-skips when the TLE data is absent.
//
// Tag deliberately avoids the substring "tle" so it is not caught by the CI
// exclude regex.

namespace {
  // A well-formed pair of TLE data lines. The first argument to FromLines is the
  // object name line; the next two are the standard TLE line 1 and line 2.
  const std::string kName = "GPS BIIR-2  (PRN 13)";
  const std::string kLine1
      = "1 25544U 98067A   24001.50000000  .00001000  00000-0  20000-3 0  9990";
  const std::string kLine2
      = "2 25544  51.6400 208.0000 0006000  90.0000 270.0000 15.50000000123456";
}  // namespace

TEST_CASE("interfaces.two_line_element.from_lines_decodes_fixed_columns") {
  TLE tle = TLE::FromLines(kName, kLine1, kLine2);

  // Epoch: two-digit year 24 -> 2024 (parse_year threshold), day-of-year 1.5.
  REQUIRE(tle.epoch_year == 2024);
  REQUIRE_THAT(tle.epoch_day, WithinAbs(1.5, 1e-6));

  // Orbital elements read straight from their fixed TLE columns.
  REQUIRE_THAT(tle.inclination, WithinAbs(51.64, 1e-6));
  REQUIRE_THAT(tle.raan, WithinAbs(208.0, 1e-6));
  REQUIRE_THAT(tle.eccentricity, WithinAbs(0.0006, 1e-9));
  REQUIRE_THAT(tle.arg_perigee, WithinAbs(90.0, 1e-6));
  REQUIRE_THAT(tle.mean_anomaly, WithinAbs(270.0, 1e-6));
  REQUIRE_THAT(tle.mean_motion, WithinAbs(15.5, 1e-6));

  // The GPS name resolves to its PRN, and the epoch converts to a finite TAI.
  REQUIRE(tle.prn == 13);
  REQUIRE(tle.name.find("GPS") == 0);
  REQUIRE(std::isfinite(tle.epoch_tai));
  REQUIRE(std::isfinite(tle.bstar));
}

TEST_CASE("interfaces.two_line_element.get_prn_branches") {
  // Names shorter than three characters have no recognizable prefix.
  REQUIRE(GetPrn("AB") == -1);

  // The first-flight Galileo (…01…) is explicitly rejected.
  REQUIRE(GetPrn("GSAT0101 (GALILEO 1)") == -1);

  // A BeiDou name whose parenthesized token does not start with 'C' maps to 0.
  REQUIRE(GetPrn("BEIDOU-3 IGSO-1 (55)") == 0);

  // QZSS PRN is taken from the digit at offset 4 of the name.
  REQUIRE(GetPrn("QZS-3 (QZSS-3)") == 3);

  // COSMOS (GLONASS) slot number is the first three digits after '('.
  REQUIRE(GetPrn("COSMOS 2501 (702K)") == 702);
}

TEST_CASE("interfaces.two_line_element.get_prn_unknown_prefix_throws") {
  // A long-enough name with an unrecognized prefix hits the terminal check.
  REQUIRE_THROWS(GetPrn("FOOBAR SATELLITE 1"));
}

TEST_CASE("interfaces.two_line_element.from_file_and_loader") {
  // Bundled GNSS TLE snapshot; skip cleanly when the data is not present.
  std::filesystem::path tle_path = GetDataPath() / "tle" / "2023_06_09" / "gps_2023_06_09.txt";
  if (!std::filesystem::exists(tle_path)) {
    SKIP("bundled TLE data not available");
  }

  const std::string basename = "gps_2023_06_09.txt";
  std::vector<TLE> tles = TLE::FromFile(basename);
  REQUIRE_FALSE(tles.empty());
  for (const TLE& t : tles) {
    REQUIRE(t.prn > 0);  // only valid PRNs are kept
    REQUIRE(t.inclination >= 0.0);
    REQUIRE(t.inclination <= 180.0);
    REQUIRE(t.eccentricity >= 0.0);
    REQUIRE(t.eccentricity < 1.0);
    REQUIRE(t.mean_motion > 0.0);
    REQUIRE(std::isfinite(t.epoch_tai));
  }
}
