#include <lupnt/dynamics/analytical_orbit_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.analytical_orbit_dynamics") {
  ClassicalOE coe(Vec6(7000.0e3, 0.01, 0.2, 0.3, 0.4, 0.5), Frame::GCRF);
  KeplerianDynamics<ClassicalOE> dyn(GM_EARTH);
  MatXd stm;
  State xf = dyn.Propagate(coe, 0.0, 60.0, nullptr, &stm);

  Real n = sqrt(GM_EARTH / pow(coe.a(), 3));
  REQUIRE_THAT(xf(0).val(), WithinAbs(coe(0).val(), epsilon));
  REQUIRE_THAT(xf(5).val(), WithinAbs(WrapToPi(coe.M() + n * 60.0).val(), epsilon));
  REQUIRE(stm.rows() == 6);
  REQUIRE(stm.cols() == 6);
  REQUIRE_THAT(stm(5, 0), WithinAbs((-1.5 * (n / coe.a() * 60.0)).val(), 1e-16));

  VecX times(3);
  times << 0.0, 10.0, 20.0;
  MatX history = dyn.Propagate(coe, times);
  REQUIRE(history.rows() == 3);
  REQUIRE(history.cols() == 6);
  REQUIRE_THAT(history(0, 5).val(), WithinAbs(coe.M().val(), epsilon));
}
