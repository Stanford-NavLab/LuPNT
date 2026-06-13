#include <lupnt/dynamics/imu_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.imu_dynamics") {
  SECTION("configuration and setters select IMU model") {
    Config config = YAML::Load("{model: LN200S}");
    ImuDynamics dyn(config);

    REQUIRE(dyn.GetModel() == ImuModel::LN200S);

    dyn.SetModel(ImuModel::UNDEFINED);
    REQUIRE(dyn.GetModel() == ImuModel::UNDEFINED);
  }

  SECTION("propagation requires a supported model") {
    ImuDynamics dyn;
    ImuState x0;

    REQUIRE_THROWS(dyn.Propagate(x0, 0.0, 1.0));
  }
}
