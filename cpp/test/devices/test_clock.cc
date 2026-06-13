#include <lupnt/agents/agent.h>
#include <lupnt/core/constants.h>
#include <lupnt/devices/clock.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

class ConstantStateAgent : public Agent {
private:
  Cart6 rv_;

public:
  explicit ConstantStateAgent(const Cart6& rv) : rv_(rv) {}

  Cart6 GetStateAt(Real t) const override {
    (void)t;
    return rv_;
  }
};

TEST_CASE("devices.clock") {
  SECTION("default Clock owns dynamics and reads a static bias") {
    Clock clock;
    ClockState3 state;
    state.b() = 1.5;
    clock.SetState(state);

    REQUIRE(clock.GetDynamics() != nullptr);
    REQUIRE(clock.GetModel() == ClockModel::UNDEFINED);
    REQUIRE_THAT(clock.Read(10.0).val(), WithinAbs(11.5, epsilon));
    REQUIRE_THAT(clock.GetTime().val(), WithinAbs(10.0, epsilon));
  }

  SECTION("config constructor sets clock model") {
    Config config = YAML::Load("{name: clk, model: USO}");
    Clock clock(config);

    REQUIRE(clock.GetName() == "clk");
    REQUIRE(clock.GetModel() == ClockModel::USO);
    REQUIRE(clock.GetDynamics() != nullptr);
    REQUIRE(clock.GetDynamics()->GetModel() == ClockModel::USO);
  }

  SECTION("clock bias can be stored in range-equivalent meters") {
    Clock clock;
    clock.SetClockBiasUnit(ClockBiasUnit::METERS);
    ClockState3 state;
    state.b() = C * 1.5;
    clock.SetState(state);

    REQUIRE(clock.GetClockBiasUnit() == ClockBiasUnit::METERS);
    REQUIRE(clock.GetState().GetUnits() == std::vector<std::string>{"m", "m/s", "m/s^2"});
    REQUIRE_THAT(clock.Read(10.0).val(), WithinAbs(11.5, epsilon));
  }

  SECTION("config constructor sets clock bias units") {
    Config config = YAML::Load("{name: clk, model: USO, clock_bias_unit: KILOMETERS}");
    Clock clock(config);

    REQUIRE(clock.GetClockBiasUnit() == ClockBiasUnit::KILOMETERS);
    REQUIRE(clock.GetDynamics()->GetClockBiasUnit() == ClockBiasUnit::KILOMETERS);
    REQUIRE(clock.GetState().GetUnits() == std::vector<std::string>{"km", "km/s", "km/s^2"});
  }

  SECTION("relativistic clock reading queries the owning agent state") {
    Vec6 rv_vec;
    rv_vec << R_EARTH, 0.0, 0.0, 0.0, 0.0, 0.0;
    ConstantStateAgent agent(Cart6(rv_vec, Frame::GCRF));
    Ptr<Clock> clock = MakePtr<Clock>();
    clock->SetUseRelativity(true);
    clock->SetRelativityCenterBody(BodyId::EARTH);
    agent.AddDevice(clock);

    Real expected_rate
        = ClockDynamics::RelativisticRateCorrection(rv_vec.head<3>(), rv_vec.tail<3>(), GM_EARTH);
    Real reading = clock->Read(10.0);

    REQUIRE(clock->GetUseRelativity());
    REQUIRE_THAT(reading.val(), WithinAbs((10.0 + expected_rate * 10.0).val(), 1.0e-15));
  }

  SECTION("relativistic clock reading supports range-equivalent bias states") {
    Vec6 rv_vec;
    rv_vec << R_EARTH, 0.0, 0.0, 0.0, 0.0, 0.0;
    ConstantStateAgent agent(Cart6(rv_vec, Frame::GCRF));
    Ptr<Clock> clock = MakePtr<Clock>();
    clock->SetClockBiasUnit(ClockBiasUnit::METERS);
    clock->SetUseRelativity(true);
    clock->SetRelativityCenterBody(BodyId::EARTH);
    agent.AddDevice(clock);

    Real expected_rate
        = ClockDynamics::RelativisticRateCorrection(rv_vec.head<3>(), rv_vec.tail<3>(), GM_EARTH);
    Real reading = clock->Read(10.0);

    REQUIRE_THAT(reading.val(), WithinAbs((10.0 + expected_rate * 10.0).val(), 1.0e-15));
    REQUIRE_THAT(clock->GetState().b().val(), WithinAbs((C * expected_rate * 10.0).val(), 1.0e-6));
  }
}
