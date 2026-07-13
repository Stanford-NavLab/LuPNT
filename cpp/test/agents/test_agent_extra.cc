#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <stdexcept>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Minimal concrete Agent: the base class has a pure-virtual GetStateAt.
  class BareAgent : public Agent {
  public:
    Cart6 GetStateAt(Real t) const override {
      (void)t;
      return Cart6(Vec6::Zero(), Frame::MOON_CI);
    }
  };

  // Minimal concrete Application: Step is pure virtual.
  class TestApp : public Application {
  public:
    void Step(Real t) override { (void)t; }
  };

  Ptr<TestApp> MakeApp(const std::string& name) {
    auto app = MakePtr<TestApp>();
    app->SetName(name);
    return app;
  }
}  // namespace

TEST_CASE("agents.agent_extra.base_accessors") {
  BareAgent agent;

  SECTION("name / time round-trip through the base setters") {
    agent.SetName("bare");
    REQUIRE(agent.GetName() == "bare");
    agent.SetTime(12.5);
    REQUIRE(agent.GetTime() == Real(12.5));
  }

  SECTION("base Agent has no intrinsic state (always empty)") {
    // Base Agent::GetState() returns an empty State() regardless of SetState.
    State rv = Cart6(Vec6(1, 2, 3, 4, 5, 6), Frame::MOON_CI);
    agent.SetState(rv);
    REQUIRE(agent.GetState().size() == 0);
  }

  SECTION("GetWorld returns nullptr and GetSimulation throws when no simulation set") {
    REQUIRE(agent.GetWorld() == nullptr);
    REQUIRE_THROWS_AS(agent.GetSimulation(), std::runtime_error);
  }
}

TEST_CASE("agents.agent_extra.devices") {
  BareAgent agent;
  agent.SetName("rover");

  // Devices are registered under device->GetName(); mimic the `<agent>/<name>`
  // convention used by the Agent(Config&) constructor.
  auto radio = MakePtr<Transmitter>();
  radio->SetName("rover/radio");
  agent.AddDevice(radio);

  SECTION("lookup by prefixed and bare name both resolve") {
    REQUIRE(agent.GetDevice("rover/radio") == radio);  // direct hit
    REQUIRE(agent.GetDevice("radio") == radio);        // via <agent>/<name>
    REQUIRE(agent.GetDevices().size() == 1);
  }

  SECTION("unknown device name throws") {
    REQUIRE_THROWS_AS(agent.GetDevice("missing"), std::runtime_error);
  }

  SECTION("adding a duplicate device name throws") {
    auto dup = MakePtr<Receiver>();
    dup->SetName("rover/radio");
    REQUIRE_THROWS_AS(agent.AddDevice(dup), std::runtime_error);
  }

  SECTION("AddDevice wires the device's back-pointer to the agent") {
    REQUIRE(radio->GetAgent() == &agent);
  }
}

TEST_CASE("agents.agent_extra.applications") {
  BareAgent agent;
  agent.SetName("lander");

  SECTION("no applications by default") {
    REQUIRE(agent.GetApplication() == nullptr);
    REQUIRE(agent.GetApplications().empty());
    REQUIRE(agent.GetApplicationByName("anything") == nullptr);
  }

  SECTION("AddApplication appends; the first is the primary; back-pointer wired") {
    auto gnc = MakeApp("lander/LanderGncApp");
    auto nav = MakeApp("lander/LanderNavApp");
    agent.AddApplication(gnc);
    agent.AddApplication(nav);

    REQUIRE(agent.GetApplications().size() == 2);
    REQUIRE(agent.GetApplication() == gnc);  // first added is primary
    REQUIRE(gnc->GetAgent() == &agent);
    REQUIRE(nav->GetAgent() == &agent);

    // Suffix match: the stored name is agent-prefixed.
    REQUIRE(agent.GetApplicationByName("LanderNavApp") == nav);
    REQUIRE(agent.GetApplicationByName("LanderGncApp") == gnc);
    // Full-name exact match also works.
    REQUIRE(agent.GetApplicationByName("lander/LanderGncApp") == gnc);
    // No match -> nullptr.
    REQUIRE(agent.GetApplicationByName("SomeOtherApp") == nullptr);
  }

  SECTION("AddApplication ignores a null application") {
    agent.AddApplication(nullptr);
    REQUIRE(agent.GetApplications().empty());
    REQUIRE(agent.GetApplication() == nullptr);
  }

  SECTION("SetApplication replaces any existing applications") {
    agent.AddApplication(MakeApp("lander/A"));
    agent.AddApplication(MakeApp("lander/B"));
    REQUIRE(agent.GetApplications().size() == 2);

    auto sole = MakeApp("lander/Sole");
    agent.SetApplication(sole);
    REQUIRE(agent.GetApplications().size() == 1);
    REQUIRE(agent.GetApplication() == sole);
    REQUIRE(sole->GetAgent() == &agent);
    REQUIRE(agent.GetApplicationByName("Sole") == sole);
    REQUIRE(agent.GetApplicationByName("A") == nullptr);
  }

  SECTION("SetApplication with null clears all applications") {
    agent.AddApplication(MakeApp("lander/A"));
    agent.SetApplication(nullptr);
    REQUIRE(agent.GetApplications().empty());
    REQUIRE(agent.GetApplication() == nullptr);
  }
}
