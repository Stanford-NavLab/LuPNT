#include <lupnt/agents/ground_station.h>
#include <lupnt/interfaces/spice.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

namespace {
  spice::GroundStationSpiceData MakeDss14Data() {
    BodyData earth = GetBodyData(BodyId::EARTH);
    LatLonAlt lla(Vec3(35.426456, 243.110461, 1001.39), earth.fixed_frame);
    Cart3 position = LatLonAltToCart(lla, earth.R, earth.flattening);

    spice::GroundStationSpiceData data;
    data.name = "DSS14";
    data.body_id = BodyId::EARTH;
    data.frame = "ITRF";
    data.position_m = position.head(3).cast<double>();
    data.latitude_deg = 35.426456;
    data.longitude_deg = 243.110461;
    data.altitude_m = 1001.39;
    return data;
  }

  void RequireStatePositionNear(const Cart6& state, const Vec3& expected, double tol) {
    for (int i = 0; i < 3; ++i) {
      REQUIRE_THAT(state.r()(i).val(), WithinAbs(expected(i).val(), tol));
    }
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(state.v()(i).val(), WithinAbs(0.0, tol));
  }
}  // namespace

TEST_CASE("agents.ground_station constructs static cartesian state from SPICE station data") {
  spice::GroundStationSpiceData data = MakeDss14Data();
  GroundStation station(data);

  REQUIRE(station.GetName() == "DSS14");
  REQUIRE(station.GetBodyId() == BodyId::EARTH);
  REQUIRE(station.GetFrame() == Frame::ITRF);
  REQUIRE_THAT(station.GetLatitudeDegDouble(), WithinAbs(data.latitude_deg, epsilon));
  REQUIRE_THAT(station.GetLongitudeDegDouble(), WithinAbs(data.longitude_deg, epsilon));
  REQUIRE_THAT(station.GetAltitudeMDouble(), WithinAbs(data.altitude_m, epsilon));

  Cart6 state = station.GetStateAt(0.0);
  REQUIRE(state.GetFrame() == Frame::ITRF);
  RequireStatePositionNear(state, data.position_m.cast<Real>(), epsilon);

  Cart6 propagated = station.GetStateAt(123.0);
  REQUIRE(propagated.GetFrame() == Frame::ITRF);
  RequireStatePositionNear(propagated, data.position_m.cast<Real>(), epsilon);
}

TEST_CASE("agents.ground_station constructs from config and factory") {
  Config config;
  config["name"] = "ManualStation";
  config["latitude_deg"] = 12.5;
  config["longitude_deg"] = -45.0;
  config["altitude_m"] = 250.0;

  GroundStation station(config);

  REQUIRE(station.GetName() == "ManualStation");
  REQUIRE(station.GetFrame() == Frame::ITRF);
  REQUIRE_THAT(station.GetLatitudeDouble(), WithinAbs(12.5, epsilon));
  REQUIRE_THAT(station.GetLongitudeDouble(), WithinAbs(-45.0, epsilon));
  REQUIRE_THAT(station.GetAltitudeDouble(), WithinAbs(250.0, epsilon));

  Config factory_config = config;
  Ptr<Agent> agent = AgentFactory::Create("GroundStation", factory_config);
  auto factory_station = std::dynamic_pointer_cast<GroundStation>(agent);
  REQUIRE(factory_station != nullptr);
  REQUIRE(factory_station->GetName() == "ManualStation");
  REQUIRE_THAT(factory_station->GetAltitudeMDouble(), WithinAbs(250.0, epsilon));
}

TEST_CASE("interfaces.spice resolves NAIF bodies for ground station kernels") {
  REQUIRE(spice::HasNaifBody("EARTH"));
  REQUIRE(spice::GetNaifId("EARTH") == static_cast<int>(BodyId::EARTH));
  REQUIRE(spice::GetNaifName(static_cast<int>(BodyId::EARTH)) == "EARTH");
  REQUIRE_FALSE(spice::HasNaifBody("__LUPNT_UNKNOWN_STATION__"));
}
