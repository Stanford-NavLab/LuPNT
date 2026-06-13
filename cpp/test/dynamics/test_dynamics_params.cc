#include <lupnt/dynamics/dynamics_params.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.dynamics_params") {
  SECTION("parameter-estimation enum values remain distinct") {
    REQUIRE(ParamsEstOption::TrueFixed != ParamsEstOption::Estimated);
    REQUIRE(ParamsEstOption::Estimated != ParamsEstOption::Consider);
    REQUIRE(static_cast<int>(ParamsEstOption::TrueFixed) == 0);
  }
}
