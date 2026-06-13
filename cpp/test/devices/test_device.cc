#include <lupnt/devices/device.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("devices.device") {
  SECTION("default Device stores name, agent pointer, and frequency") {
    Device device;

    device.SetName("payload");
    device.SetFrequency(2.5);
    device.SetAgent(nullptr);

    REQUIRE(device.GetName() == "payload");
    REQUIRE_THAT(device.GetFrequency().val(), WithinAbs(2.5, epsilon));
    REQUIRE(device.GetAgent() == nullptr);
  }

  SECTION("config constructor reads name and frequency") {
    Config config = YAML::Load("{name: sensor, frequency: 4.0}");
    Device device(config);

    REQUIRE(device.GetName() == "sensor");
    REQUIRE_THAT(device.GetFrequency().val(), WithinAbs(4.0, epsilon));
  }
}
