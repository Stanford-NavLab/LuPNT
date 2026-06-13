#include <lupnt/data/iau_sofa.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("data.iau_sofa") {
  IauSofaData j2000 = GetIauSofaData(JD_J2000_TT);
  REQUIRE(std::isfinite(j2000.X.val()));
  REQUIRE(std::isfinite(j2000.Y.val()));
  REQUIRE(std::isfinite(j2000.s.val()));
  REQUIRE(std::abs(j2000.X.val()) < 1.0e4);
  REQUIRE(std::abs(j2000.Y.val()) < 1.0e4);
  REQUIRE(std::abs(j2000.s.val()) < 1.0);
}
