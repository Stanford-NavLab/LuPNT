#include <lupnt/interfaces/spice_cheby.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("interfaces.spice_cheby") {
  double scale[2] = {5.0, 5.0};
  double coeff[2] = {10.0, 2.0};
  double f = 0.0;
  double df = 0.0;
  cheby_eval(7.5, scale, coeff, 2, &f, &df);
  REQUIRE_THAT(f, Catch::Matchers::WithinAbs(11.0, epsilon));
  REQUIRE_THAT(df, Catch::Matchers::WithinAbs(0.4, epsilon));

  Vec2 ad_value = cheby_eval_ad(7.5, scale, coeff, 2);
  REQUIRE_THAT(ad_value(0).val(), Catch::Matchers::WithinAbs(11.0, epsilon));
  REQUIRE_THAT(ad_value(1).val(), Catch::Matchers::WithinAbs(0.4, epsilon));

  double seg[12] = {5.0, 5.0, 10.0, 2.0, 20.0, -1.0, 30.0, 0.5, 0.0, 10.0, 8.0, 1.0};
  REQUIRE(cheby_verify(seg, 12) == 0);
  double pos[3] = {};
  double vel[3] = {};
  REQUIRE(cheby_posvel(7.5, seg, 12, pos, vel) == 0);
  REQUIRE_THAT(pos[0], Catch::Matchers::WithinAbs(11.0, epsilon));
  REQUIRE_THAT(pos[1], Catch::Matchers::WithinAbs(19.5, epsilon));
  REQUIRE_THAT(pos[2], Catch::Matchers::WithinAbs(30.25, epsilon));
  REQUIRE_THAT(vel[0], Catch::Matchers::WithinAbs(0.4, epsilon));
  REQUIRE_THAT(vel[1], Catch::Matchers::WithinAbs(-0.2, epsilon));
  REQUIRE_THAT(vel[2], Catch::Matchers::WithinAbs(0.1, epsilon));
  REQUIRE(cheby_posvel(-1.0, seg, 12, pos, vel) == 1);
}
