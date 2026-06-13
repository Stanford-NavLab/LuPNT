#include <lupnt/devices/imu.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("devices.imu") {
  SECTION("LN200S parameters are finite and positive") {
    auto params = Imu::GetParameters(ImuModel::LN200S);

    REQUIRE(params.sigma_a.val() > 0.0);
    REQUIRE(params.sigma_a_bias.val() > 0.0);
    REQUIRE(params.sigma_w.val() > 0.0);
    REQUIRE(params.sigma_w_bias.val() > 0.0);
  }

  SECTION("data buffer can be emptied") {
    Imu imu;

    REQUIRE(imu.GetData().empty());
    imu.EmptyData();
    REQUIRE(imu.GetData().empty());
  }
}
