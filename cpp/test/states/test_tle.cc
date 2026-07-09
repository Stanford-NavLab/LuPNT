#include <lupnt/interfaces/tle.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("states.tle") {
  REQUIRE(GetPrn("GPS BIIR-2  (PRN 13)") == 13);
  REQUIRE(GetPrn("BEIDOU-2 M3 (C11)") == 11);
  REQUIRE(GetPrn("GSAT0202 (GALILEO 6)") == 6);
  REQUIRE(GetPrn("COSMOS 2501 (702K)") == 702);
  REQUIRE(GetPrn("QZS-1R (QZSS/PRN 196)") == 1);

  std::ifstream file(GetFilePath("gps_2023_06_09.txt"));
  std::string name, line2, line3;
  std::getline(file, name);
  std::getline(file, line2);
  std::getline(file, line3);

  TLE tle = TLE::FromLines(name, line2, line3);
  REQUIRE(tle.name.find("GPS") == 0);
  REQUIRE(tle.prn == 13);
  REQUIRE(tle.epoch_year == 2023);
  REQUIRE(tle.eccentricity > 0.0);
  REQUIRE(tle.mean_motion > 0.0);
}
