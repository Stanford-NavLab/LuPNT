#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

namespace {
  class TestAgent : public AgentWithDynamics {
  public:
    Cart6 GetStateAt(Real t) const override {
      (void)t;
      return GetState();
    }
  };
}  // namespace

TEST_CASE("agents.agent") {
  TestAgent agent;

  // Time
  agent.SetTime(1.0);
  REQUIRE(agent.GetTime() == 1.0);

  // Name
  agent.SetName("agent");
  REQUIRE(agent.GetName() == "agent");

  // State
  State coe = ClassicalOE({7000.0, 0.001, 0.0, 0.0, 0.0, 0.0}, Frame::MOON_CI);
  State rv = ClassicalToCart(coe, GM_MOON);
  agent.SetState(rv);
  REQUIRE(agent.GetState() == rv);

  // Attitude
  Attitude att(Vec4(1.0, 0.0, 0.0, 0.0), Vec3::Zero(), Frame::MOON_CI);
  agent.SetAttitude(att);
  REQUIRE(agent.GetAttitude() == att);

  // Dynamics
  auto dyn_nbody = std::make_shared<NBodyDynamics>();
  dyn_nbody->AddBody(Body::Moon());
  dyn_nbody->SetFrame(Frame::MOON_CI);

  Ptr<Dynamics> dyn = dyn_nbody;
  agent.SetDynamics(dyn);
  REQUIRE(agent.GetDynamics() == dyn.get());

  // Devices
  auto tx = MakePtr<Transmitter>();
  tx->SetName("tx");
  auto rx = MakePtr<Receiver>();
  rx->SetName("rx");
  auto transponder = MakePtr<Transponder>();
  transponder->SetName("transponder");
  agent.AddDevice(tx);
  agent.AddDevice(rx);
  agent.AddDevice(transponder);
  REQUIRE(agent.GetDevices().size() == 3);
  REQUIRE(agent.GetDevice("tx") == tx);
  REQUIRE(agent.GetDevice("rx") == rx);
  REQUIRE(agent.GetDevice("transponder") == transponder);
}
