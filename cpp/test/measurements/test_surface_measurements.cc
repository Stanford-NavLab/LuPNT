#include <lupnt/core/constants.h>
#include <lupnt/measurements/surface_measurements.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("measurements.surface_lans") {
  SurfaceLansMeasurement m;
  m.r_sat = Vec3d(1000.0, 0.0, 0.0);
  m.sigma_m = 1.0;
  m.sise_m = 3.0;

  Vec3d r_rover(0.0, 0.0, 0.0);
  const double cb = 2.0e-6;

  REQUIRE_THAT(m.PredictedRange(r_rover, cb), WithinRel(1000.0 + C * cb, 1e-12));
  REQUIRE_THAT(m.NoiseVariance(), WithinAbs(1.0 + 9.0, 1e-12));

  // Line-of-sight unit vector points from the rover toward the satellite (+x here).
  Vec3d u = m.LosUnit(r_rover);
  REQUIRE_THAT(u(0), WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(u(1), WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(u(2), WithinAbs(0.0, 1e-12));

  // The range partial d(rho)/d(r_rover) equals -u; verify against a finite difference.
  const double eps = 1e-3;
  for (int i = 0; i < 3; ++i) {
    Vec3d rp = r_rover, rm = r_rover;
    rp(i) += eps;
    rm(i) -= eps;
    double fd = (m.PredictedRange(rp, cb) - m.PredictedRange(rm, cb)) / (2 * eps);
    REQUIRE_THAT(-u(i), WithinAbs(fd, 1e-5));
  }
}
