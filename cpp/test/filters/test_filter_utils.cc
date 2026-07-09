#include <lupnt/numerics/filters/filter_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("filters.filter_utils") {
  MatXd P = InitialCovariancePosVelClock(2.0, 3.0, 4.0, 5.0);
  REQUIRE(P.rows() == 8);
  REQUIRE(P.cols() == 8);
  REQUIRE_THAT(P(0, 0), WithinAbs(4.0, epsilon));
  REQUIRE_THAT(P(3, 3), WithinAbs(9.0, epsilon));
  REQUIRE_THAT(P(6, 6), WithinAbs(16.0, epsilon));
  REQUIRE_THAT(P(7, 7), WithinAbs(25.0, epsilon));

  Vec3d coeffs = ProcessNoisePosVelCoeffs(2.0);
  REQUIRE_THAT(coeffs(0), WithinAbs(8.0 / 3.0, epsilon));
  REQUIRE_THAT(coeffs(1), WithinAbs(2.0, epsilon));
  REQUIRE_THAT(coeffs(2), WithinAbs(2.0, epsilon));

  MatXd Q = ProcessNoisePosVel(Mat3d::Identity(), 2.0);
  REQUIRE(Q.rows() == 6);
  REQUIRE(Q.cols() == 6);
  REQUIRE_THAT(Q(0, 3), WithinAbs(2.0, epsilon));

  MatXd Phi = StateTransitionMatrixPosVel(5.0, 3);
  REQUIRE_THAT(Phi(0, 3), WithinAbs(5.0, epsilon));
  REQUIRE_THAT(Phi(5, 5), WithinAbs(1.0, epsilon));
}
