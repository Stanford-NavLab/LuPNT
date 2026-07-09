#include <lupnt/simulations/simulation.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.simulation") {
  SECTION("scheduled events run in time order through the configured duration") {
    Simulation sim;
    sim.SetDuration(2.0);

    std::vector<double> calls;
    sim.Schedule(1.0, [&](Real t) { calls.push_back(t.val()); });
    sim.Schedule(0.5, [&](Real t) { calls.push_back(t.val()); });
    sim.Schedule(3.0, [&](Real t) { calls.push_back(t.val()); });

    sim.Run();

    REQUIRE(calls == std::vector<double>{0.5, 1.0});
    REQUIRE_THAT(sim.GetTime().val(), WithinAbs(1.0, epsilon));
  }

  SECTION("Publish schedules subscriber callbacks") {
    Simulation sim;
    sim.SetDuration(1.0);

    std::string received;
    sim.Subscribe("topic",
                  [&](const std::any& payload) { received = std::any_cast<std::string>(payload); });
    sim.Publish(0.25, "topic", std::string("payload"));

    sim.Run();

    REQUIRE(received == "payload");
  }
}
