#include <lupnt/devices/camera.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("devices.camera") {
  SECTION("Camera can be constructed from config and stepped") {
    Config config = YAML::Load("{name: camera, frequency: 1.0}");
    Camera camera(config);

    REQUIRE(camera.GetName() == "camera");
    REQUIRE_THAT(camera.GetFrequency().val(), WithinAbs(1.0, epsilon));
    REQUIRE_NOTHROW(camera.Setup());
    REQUIRE_NOTHROW(camera.Step(0.0));
  }
}
