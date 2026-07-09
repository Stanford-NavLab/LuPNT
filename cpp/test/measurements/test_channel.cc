#include <lupnt/devices/comm_devices.h>
#include <lupnt/measurements/channel.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

namespace {
  class TestReceiver : public Receiver {
  public:
    void Receive(Real t, void* data) override {
      last_t = t;
      last_data = data;
      ++count;
    }

    int count = 0;
    Real last_t = 0.0;
    void* last_data = nullptr;
  };
}  // namespace

TEST_CASE("measurements.channel") {
  Config config = YAML::Load("name: unit-test-channel\nclass: Channel\n");
  Channel channel(config);
  REQUIRE(channel.GetName() == "unit-test-channel");

  Transmitter tx;
  TestReceiver rx;
  channel.AddDevice(&tx);
  channel.AddDevice(&rx);

  int payload = 42;
  channel.Send(&tx, 12.5, &payload);

  REQUIRE(rx.count == 1);
  REQUIRE_THAT(rx.last_t.val(), Catch::Matchers::WithinAbs(12.5, epsilon));
  REQUIRE(rx.last_data == &payload);
  REQUIRE_THROWS(channel.Send(&tx, 0.0, nullptr));

  auto from_config = Channel::FromConfig(config);
  REQUIRE(from_config->GetName() == "unit-test-channel");
}
