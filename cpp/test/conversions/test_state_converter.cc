#include <lupnt/conversions/state_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.state_converter") {
  SECTION("ConvertState dispatches absolute conversions by type") {
    ClassicalOE coe(Vec6(7000.0e3, 0.01, 20.0 * RAD, 30.0 * RAD, 40.0 * RAD, 50.0 * RAD),
                    Frame::GCRF);

    State rv = ConvertState(coe, Cart6::TYPE, GM_EARTH);
    State recovered = ConvertState(rv, ClassicalOE::TYPE, GM_EARTH);

    REQUIRE(rv.GetType() == Cart6::TYPE);
    REQUIRE(recovered.GetType() == ClassicalOE::TYPE);
    REQUIRE_THAT(recovered(0).val(), WithinAbs(coe(0).val(), 1.0e-6));
    REQUIRE_THAT(recovered(1).val(), WithinAbs(coe(1).val(), 1.0e-12));
  }

  SECTION("requesting the same state type returns an equivalent state") {
    Cart6 rv(Vec6(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), Frame::GCRF);

    State out = ConvertState(rv, Cart6::TYPE);

    REQUIRE(out.GetType() == Cart6::TYPE);
    REQUIRE(out.GetFrame() == Frame::GCRF);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(out(i).val(), WithinAbs(rv(i).val(), epsilon));
  }
}
