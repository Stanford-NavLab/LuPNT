#include <lupnt/conversions/coordinate_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.coordinate_conversions") {
  SECTION("latitude longitude altitude converts to Cartesian and back") {
    LatLonAlt lla(Vec3(35.0 * RAD, -120.0 * RAD, 250.0), Frame::ITRF);

    State xyz = LatLonAltToCart(lla, R_EARTH, WGS84_F);
    State recovered = CartToLatLonAlt(xyz, R_EARTH, WGS84_F);

    REQUIRE(recovered.GetType() == LatLonAlt::TYPE);
    REQUIRE_THAT(recovered(0).val(), WithinAbs(lla(0).val(), 1.0e-10));
    REQUIRE_THAT(recovered(1).val(), WithinAbs(lla(1).val(), 1.0e-10));
    REQUIRE_THAT(recovered(2).val(), WithinAbs(lla(2).val(), 1.0e-5));
  }

  SECTION("azimuth elevation range converts to ENU and back") {
    AzElRange aer(Vec3(30.0 * RAD, 20.0 * RAD, 1234.5));

    State enu = AzElRangeToEastNorthUp(aer);
    State recovered = EastNorthUpToAzElRange(enu);

    REQUIRE(recovered.GetType() == AzElRange::TYPE);
    REQUIRE_THAT(recovered(0).val(), WithinAbs(aer(0).val(), epsilon));
    REQUIRE_THAT(recovered(1).val(), WithinAbs(aer(1).val(), epsilon));
    REQUIRE_THAT(recovered(2).val(), WithinAbs(aer(2).val(), epsilon));
  }

  SECTION("ENU to Cartesian round trip around a reference point") {
    LatLonAlt lla(Vec3(10.0 * RAD, 20.0 * RAD, 0.0), Frame::ITRF);
    State ref = LatLonAltToCart(lla, R_EARTH, WGS84_F);
    Cart3 enu(Vec3(100.0, 200.0, 300.0), Frame::ITRF);

    State xyz = EastNorthUpToCart(enu, ref, R_EARTH, WGS84_F);
    State recovered = CartToEastNorthUp(xyz, ref, R_EARTH, WGS84_F);

    for (int i = 0; i < 3; ++i) REQUIRE_THAT(recovered(i).val(), WithinAbs(enu(i).val(), 1.0e-6));
  }
}
