#include <lupnt/conversions/state_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.state_conversions") {
  SECTION("Classical orbital elements convert to Cartesian and back") {
    ClassicalOE coe(Vec6(7000.0e3, 0.05, 40.0 * RAD, 30.0 * RAD, 20.0 * RAD, 10.0 * RAD),
                    Frame::GCRF);

    State rv = ClassicalToCart(coe, GM_EARTH);
    State recovered = CartToClassical(rv, GM_EARTH);

    REQUIRE(rv.GetType() == Cart6::TYPE);
    REQUIRE(rv.GetFrame() == coe.GetFrame());
    REQUIRE(recovered.GetType() == ClassicalOE::TYPE);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(recovered(i).val(), WithinAbs(coe(i).val(), 1.0e-6));
  }

  SECTION("alternate orbital element sets round trip through ClassicalOE") {
    ClassicalOE coe(Vec6(8000.0e3, 0.1, 35.0 * RAD, 40.0 * RAD, 50.0 * RAD, 60.0 * RAD),
                    Frame::GCRF);

    State qn = ClassicalToQuasiNonsing(coe, GM_EARTH);
    State eq = ClassicalToEquinoctial(coe, GM_EARTH);
    State del = ClassicalToDelaunay(coe, GM_EARTH);

    State coe_from_qn = QuasiNonsingToClassical(qn, GM_EARTH);
    State coe_from_eq = EquinoctialToClassical(eq, GM_EARTH);
    State coe_from_del = DelaunayToClassical(del, GM_EARTH);

    REQUIRE(qn.GetType() == QuasiNonsingularOE::TYPE);
    REQUIRE(eq.GetType() == EquinoctialOE::TYPE);
    REQUIRE(del.GetType() == DelaunayOE::TYPE);
    REQUIRE_THAT(coe_from_qn(0).val(), WithinAbs(coe(0).val(), 1.0e-6));
    REQUIRE_THAT(coe_from_eq(1).val(), WithinAbs(coe(1).val(), 1.0e-12));
    REQUIRE_THAT(coe_from_del(2).val(), WithinAbs(coe(2).val(), 1.0e-12));
  }
}
