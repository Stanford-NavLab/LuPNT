#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "lupnt/core/file.h"
#include "lupnt/environment/plasma/core/string_file_utils.h"
#include "lupnt/environment/plasma/core/user_filepath.h"

using namespace pecsim;
using namespace Catch::Matchers;

// The plasma module locates its runtime data relative to a base path. In the
// pixi test environment PECSIMPY_BASE_PATH is exported already, but we set it
// explicitly here so the test is self-contained.
static void SetupPlasmaBasePath() {
  std::filesystem::path plasma_base = lupnt::GetDataPath() / "plasma";
  set_base_path(plasma_base.string());
}

// ------------------------------ trim_right_copy ----------------------------

TEST_CASE("environment.plasma.string_file_utils.trim_right_copy") {
  SECTION("trailing spaces are removed") { REQUIRE(trim_right_copy("hello   ") == "hello"); }
  SECTION("mixed trailing whitespace is removed") {
    REQUIRE(trim_right_copy("value\t \n") == "value");
  }
  SECTION("leading and internal whitespace is preserved") {
    REQUIRE(trim_right_copy("  a b  ") == "  a b");
  }
  SECTION("string with no trailing whitespace is unchanged") {
    REQUIRE(trim_right_copy("abc") == "abc");
  }
  SECTION("empty and all-whitespace strings collapse to empty") {
    REQUIRE(trim_right_copy("").empty());
    REQUIRE(trim_right_copy("    ").empty());
  }
}

// ------------------------------- split_string ------------------------------

TEST_CASE("environment.plasma.string_file_utils.split_string") {
  SECTION("splits comma-separated tokens") {
    std::vector<std::string> parts = split_string("a,b,c", ',');
    REQUIRE(parts.size() == 3);
    REQUIRE(parts[0] == "a");
    REQUIRE(parts[1] == "b");
    REQUIRE(parts[2] == "c");
  }
  SECTION("a trailing separator yields a trailing empty token") {
    std::vector<std::string> parts = split_string("a,", ',');
    REQUIRE(parts.size() == 2);
    REQUIRE(parts[0] == "a");
    REQUIRE(parts[1].empty());
  }
  SECTION("no separator returns the whole string as one token") {
    std::vector<std::string> parts = split_string("no-sep-here", ',');
    REQUIRE(parts.size() == 1);
    REQUIRE(parts[0] == "no-sep-here");
  }
  SECTION("empty input returns a single empty token") {
    std::vector<std::string> parts = split_string("", ',');
    REQUIRE(parts.size() == 1);
    REQUIRE(parts[0].empty());
  }
  SECTION("splitting on a different separator works") {
    std::vector<std::string> parts = split_string("2024 100 12", ' ');
    REQUIRE(parts.size() == 3);
    REQUIRE(parts[0] == "2024");
    REQUIRE(parts[2] == "12");
  }
}

// ---------------------------- find_file_in_dir -----------------------------

TEST_CASE("environment.plasma.string_file_utils.find_file_in_dir") {
  SetupPlasmaBasePath();
  std::filesystem::path data_dir = std::filesystem::path(get_base_path()) / "data";

  SECTION("locates a bundled IRI data file by full name") {
    std::optional<std::filesystem::path> found = find_file_in_dir(data_dir, "apf107.dat");
    REQUIRE(found.has_value());
    REQUIRE(found.value().filename().string() == "apf107.dat");
    REQUIRE(std::filesystem::exists(found.value()));
  }
  SECTION("locates a bundled IRI data file by stem") {
    std::optional<std::filesystem::path> found = find_file_in_dir(data_dir, "apf107");
    REQUIRE(found.has_value());
    REQUIRE(found.value().stem().string() == "apf107");
  }
  SECTION("returns nullopt for a file that does not exist") {
    std::optional<std::filesystem::path> found
        = find_file_in_dir(data_dir, "this_file_does_not_exist_zzz.dat");
    REQUIRE_FALSE(found.has_value());
  }
}

// NOTE: get_file_path() is intentionally not exercised here: the header
// declares get_file_path(const std::string&) but the implementation defines
// get_file_path(std::string_view), so the declared overload has no definition
// and cannot be linked from a test. find_file_in_dir() (its core logic) is
// covered above instead.

// --------------------------------- TLE parse -------------------------------

// A well-formed pair of TLE data lines with fields at the exact columns the
// parser reads. Values chosen so every parsed field is a clean, exact number.
static const std::string kTleLine2
    = "1 25544U 98067A   24001.50000000  .00001000  00000-0  20000-3 0  9990";
static const std::string kTleLine3
    = "2 25544  51.6400 208.0000 0006000  90.0000 270.0000 15.50000000123456";

TEST_CASE("environment.plasma.string_file_utils.tle_from_lines") {
  SECTION("parses orbital elements from fixed columns") {
    TLE tle = TLE::FromLines("ISS (ZARYA)", kTleLine2, kTleLine3);
    REQUIRE_THAT(tle.epoch_year, WithinAbs(24.0, 1e-9));
    REQUIRE_THAT(tle.epoch_day, WithinAbs(1.5, 1e-6));
    REQUIRE_THAT(tle.inclination, WithinAbs(51.64, 1e-6));
    REQUIRE_THAT(tle.raan, WithinAbs(208.0, 1e-6));
    REQUIRE_THAT(tle.eccentricity, WithinAbs(0.0006, 1e-9));
    REQUIRE_THAT(tle.arg_perigee, WithinAbs(90.0, 1e-6));
    REQUIRE_THAT(tle.mean_anomaly, WithinAbs(270.0, 1e-6));
    REQUIRE_THAT(tle.mean_motion, WithinAbs(15.5, 1e-6));
    // A name without a recognized GNSS prefix maps to prn == -1.
    REQUIRE(tle.prn == -1);
    // Epoch was converted to a finite TAI-seconds value.
    REQUIRE(std::isfinite(tle.epoch_utc));
  }
  SECTION("extracts the PRN from a GPS satellite name") {
    TLE tle = TLE::FromLines("GPS BIIR-2  (PRN 13)", kTleLine2, kTleLine3);
    REQUIRE(tle.prn == 13);
    // Orbital elements parse identically regardless of the name field.
    REQUIRE_THAT(tle.inclination, WithinAbs(51.64, 1e-6));
  }
}
