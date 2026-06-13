#include <lupnt/core/constants.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("core.constants.units") {
  SECTION("unit systems scale SI quantities by dimensional powers") {
    REQUIRE_THAT(KM_S_KG_UNITS.Length(R_EARTH), WithinRel(R_EARTH / 1000.0, 1.0e-12));
    REQUIRE_THAT(KM_S_KG_UNITS.Velocity(C), WithinRel(C / 1000.0, 1.0e-12));
    REQUIRE_THAT(KM_S_KG_UNITS.GravitationalParameter(GM_EARTH),
                 WithinRel(GM_EARTH / 1.0e9, 1.0e-12));
    REQUIRE_THAT(KM_S_KG_UNITS.Pressure(P_SUN), WithinRel(P_SUN * 1000.0, 1.0e-12));
  }

  SECTION("physical constants can be requested in kilometers seconds kilograms") {
    PhysicalConstants km = GetPhysicalConstants(KM_S_KG_UNITS);

    REQUIRE_THAT(km.R_EARTH, WithinRel(6378.137, 1.0e-12));
    REQUIRE_THAT(km.GM_EARTH, WithinRel(398600.435507, 1.0e-12));
    REQUIRE_THAT(km.AU, WithinRel(AU / 1000.0, 1.0e-12));
    REQUIRE_THAT(km.C, WithinRel(C / 1000.0, 1.0e-12));
  }

  SECTION("coordinate scales apply IAU constant scaling separately from units") {
    REQUIRE(AreCoordinateScalesConvertible(CoordinateScale::TCB, CoordinateScale::TDB));
    REQUIRE(AreCoordinateScalesConvertible(CoordinateScale::TCG, CoordinateScale::TT));
    REQUIRE(AreCoordinateScalesConvertible(CoordinateScale::TCL, CoordinateScale::TL));
    REQUIRE_FALSE(AreCoordinateScalesConvertible(CoordinateScale::TDB, CoordinateScale::TT));

    REQUIRE_THAT(CheckedCoordinateScaleRatio(CoordinateScale::TCB, CoordinateScale::TDB),
                 WithinRel(1.0 - L_B, 1.0e-15));
    REQUIRE_THAT(CheckedCoordinateScaleRatio(CoordinateScale::TDB, CoordinateScale::TCB),
                 WithinRel(1.0 / (1.0 - L_B), 1.0e-15));
    REQUIRE_THROWS_AS(CheckedCoordinateScaleRatio(CoordinateScale::TDB, CoordinateScale::TT),
                      std::invalid_argument);

    Vec6 rv;
    rv << 1.0, 2.0, 3.0, 0.1, 0.2, 0.3;
    Vec6 rv_tcb = ScaleStateForCoordinateScale(rv, CoordinateScale::TDB, CoordinateScale::TCB);
    REQUIRE_THAT(rv_tcb(0), WithinRel(1.0 / (1.0 - L_B), 1.0e-15));
    REQUIRE_THAT(rv_tcb(1), WithinRel(2.0 / (1.0 - L_B), 1.0e-15));
    REQUIRE_THAT(rv_tcb(3), WithinRel(0.1, 1.0e-15));
    REQUIRE_THAT(rv_tcb(5), WithinRel(0.3, 1.0e-15));
  }

  SECTION("physical constants can be requested in a supported coordinate scale") {
    PhysicalConstants tdb_km = GetPhysicalConstants(KM_S_KG_UNITS);
    PhysicalConstants tcb_km = GetPhysicalConstants(KM_S_KG_UNITS, CoordinateScale::TCB);

    REQUIRE(tdb_km.coordinate_scale == CoordinateScale::TDB);
    REQUIRE(tcb_km.coordinate_scale == CoordinateScale::TCB);
    REQUIRE_THAT(tcb_km.GM_EARTH, WithinRel(tdb_km.GM_EARTH / (1.0 - L_B), 1.0e-15));
    REQUIRE_THAT(tcb_km.AU, WithinRel(tdb_km.AU / (1.0 - L_B), 1.0e-15));
    REQUIRE_THAT(tcb_km.C, WithinRel(tdb_km.C, 1.0e-15));
    REQUIRE_THAT(tcb_km.OMEGA_EARTH, WithinRel(tdb_km.OMEGA_EARTH * (1.0 - L_B), 1.0e-15));
    REQUIRE_THROWS_AS(GetPhysicalConstants(KM_S_KG_UNITS, CoordinateScale::TCL),
                      std::invalid_argument);
  }
}
