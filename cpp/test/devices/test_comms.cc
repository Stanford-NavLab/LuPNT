#include <lupnt/devices/comm_devices.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("devices.comms") {
  SECTION("Transmitter and Receiver read config through Device") {
    Config tx_config = YAML::Load("{name: tx, frequency: 2.0}");
    Config rx_config = YAML::Load("{name: rx, frequency: 3.0}");

    Transmitter tx(tx_config);
    Receiver rx(rx_config);

    REQUIRE(tx.GetName() == "tx");
    REQUIRE(rx.GetName() == "rx");
    REQUIRE_THAT(tx.GetFrequency().val(), WithinAbs(2.0, epsilon));
    REQUIRE_THAT(rx.GetFrequency().val(), WithinAbs(3.0, epsilon));
  }

  SECTION("Receiver stores received data pointers") {
    Receiver rx;
    int payload = 7;

    REQUIRE_NOTHROW(rx.Receive(1.0, &payload));
  }
}
