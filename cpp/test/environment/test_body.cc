#include <lupnt/environment/body.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("environment.body") {
  SECTION("body data helpers expose canonical constants") {
    BodyData earth_data = GetBodyData(BodyId::EARTH);
    REQUIRE(earth_data.id == BodyId::EARTH);
    REQUIRE(earth_data.name == "EARTH");
    REQUIRE(earth_data.fixed_frame == Frame::ITRF);
    REQUIRE(earth_data.inertial_frame == Frame::GCRF);
    REQUIRE_THAT(GetBodyGM(BodyId::EARTH), WithinRel(GM_EARTH, epsilon));
    REQUIRE_THAT(GetBodyRadius(BodyId::EARTH), WithinRel(R_EARTH, epsilon));
    REQUIRE_THAT(GetBodyOmega(BodyId::EARTH), WithinRel(OMEGA_EARTH, epsilon));
    REQUIRE_THAT(GetBodyFlattening(BodyId::EARTH), WithinRel(WGS84_F, epsilon));
    REQUIRE(GetBodyName(BodyId::EARTH) == "EARTH");
    REQUIRE(GetBodyFixedFrameName(BodyId::EARTH) == Frame::ITRF);
    REQUIRE(GetInertialFrameName(BodyId::EARTH) == Frame::GCRF);
  }

  SECTION("factory constructors populate point-mass bodies without gravity fields") {
    Body moon = Body::Moon();
    REQUIRE(moon.id == BodyId::MOON);
    REQUIRE(moon.name == "MOON");
    REQUIRE(moon.fixed_frame == Frame::MOON_PA);
    REQUIRE(moon.inertial_frame == Frame::MOON_CI);
    REQUIRE_THAT(moon.GM.val(), WithinRel(GM_MOON, epsilon));
    REQUIRE_THAT(moon.R.val(), WithinRel(R_MOON, epsilon));
    REQUIRE_THAT(moon.omega.val(), WithinRel(OMEGA_MOON, epsilon));
    REQUIRE_FALSE(moon.use_gravity_field);

    Body sun = Body::Sun();
    REQUIRE(sun.id == BodyId::SUN);
    REQUIRE(sun.name == "SUN");
    REQUIRE(sun.inertial_frame == Frame::ICRF);
    REQUIRE_THAT(sun.GM.val(), WithinRel(GM_SUN, epsilon));
    REQUIRE_FALSE(sun.use_gravity_field);
  }

  SECTION("free and static CreateBody helpers agree") {
    Body mars = CreateBody(BodyId::MARS);
    Body mars_static = Body::CreateBody(BodyId::MARS);
    REQUIRE(mars.id == BodyId::MARS);
    REQUIRE(mars_static.id == BodyId::MARS);
    REQUIRE(mars.name == mars_static.name);
    REQUIRE_THAT(mars.GM.val(), WithinRel(mars_static.GM.val(), epsilon));
    REQUIRE_THAT(mars.R.val(), WithinRel(mars_static.R.val(), epsilon));
    REQUIRE_FALSE(mars.use_gravity_field);
  }

  SECTION("planet factory constants match BodyData") {
    Body venus = Body::Venus();
    BodyData venus_data = GetBodyData(BodyId::VENUS);
    REQUIRE(venus.id == venus_data.id);
    REQUIRE(venus.name == venus_data.name);
    REQUIRE(venus.fixed_frame == venus_data.fixed_frame);
    REQUIRE_THAT(venus.GM.val(), WithinRel(venus_data.GM.val(), epsilon));
    REQUIRE_THAT(venus.R.val(), WithinRel(venus_data.R.val(), epsilon));
    REQUIRE_THAT(venus.omega.val(), WithinRel(venus_data.omega.val(), epsilon));
  }

  SECTION("body helpers can expose kilometer-based constants") {
    BodyData earth_km = GetBodyData(BodyId::EARTH, KM_S_KG_UNITS);
    REQUIRE(earth_km.units == KM_S_KG_UNITS);
    REQUIRE_THAT(earth_km.GM.val(), WithinRel(GM_EARTH / 1.0e9, epsilon));
    REQUIRE_THAT(earth_km.R.val(), WithinRel(R_EARTH / 1000.0, epsilon));
    REQUIRE_THAT(GetBodyGM(BodyId::EARTH, KM_S_KG_UNITS), WithinRel(GM_EARTH / 1.0e9, epsilon));
    REQUIRE_THAT(GetBodyRadius(BodyId::EARTH, KM_S_KG_UNITS), WithinRel(R_EARTH / 1000.0, epsilon));

    Body moon_km = Body::Moon(KM_S_KG_UNITS);
    REQUIRE(moon_km.units == KM_S_KG_UNITS);
    REQUIRE_THAT(moon_km.GM.val(), WithinRel(GM_MOON / 1.0e9, epsilon));
    REQUIRE_THAT(moon_km.R.val(), WithinRel(R_MOON / 1000.0, epsilon));
  }
}
