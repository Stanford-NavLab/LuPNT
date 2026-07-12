// Unit tests for the terrain-relative lander measurement models (ex11):
// `LanderAltimeterMeasurement` (nadir height above terrain) and
// `LanderCraterMeasurement` (crater-landmark body line-of-sight). Both are
// error-state models linearized about a `NavErrorContext`; the values, noise
// covariances, and position/attitude Jacobian blocks are checked against
// hand-computed and finite-difference references. Fully self-contained.

#include <lupnt/core/constants.h>
#include <lupnt/measurements/lander_measurements.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  Mat3d SmallRotation() {
    // A modest body-to-nav rotation (~15 deg about z) to exercise the R_b2n^T terms.
    const double a = 0.26179939;  // 15 deg
    Mat3d R;
    R << std::cos(a), -std::sin(a), 0.0, std::sin(a), std::cos(a), 0.0, 0.0, 0.0, 1.0;
    return R;
  }
}  // namespace

TEST_CASE("measurements.lander_altimeter") {
  LanderAltimeterMeasurement m;
  m.altitude_m = 500.0;
  m.sigma_m = 2.0;
  m.predicted_altitude_m = 498.0;
  m.h_pos = Vec3d(0.1, -0.2, 0.97);  // DEM slope + up-axis coupling
  m.config.i_dr = 0;

  NavErrorContext nom;
  nom.r = Vec3d(10.0, 20.0, 1730000.0);
  nom.error_state_size = kLanderNavErrorStateSize;

  MatXd H;
  MeasData md = m.Compute(nom, &H);

  REQUIRE(md.value.size() == 1);
  REQUIRE_THAT(md.value(0), WithinAbs(498.0, 1e-12));  // predicted altitude
  REQUIRE_THAT(md.covariance(0, 0), WithinAbs(4.0, 1e-12));
  REQUIRE_THAT(m.NoiseVariance(), WithinAbs(4.0, 1e-12));

  REQUIRE(H.rows() == 1);
  REQUIRE(H.cols() == kLanderNavErrorStateSize);
  // Jacobian couples only to the position error block via h_pos^T.
  REQUIRE_THAT(H(0, 0), WithinAbs(0.1, 1e-12));
  REQUIRE_THAT(H(0, 1), WithinAbs(-0.2, 1e-12));
  REQUIRE_THAT(H(0, 2), WithinAbs(0.97, 1e-12));
  // Everything past the position block is zero (no velocity/attitude/bias coupling).
  REQUIRE(H.rightCols(kLanderNavErrorStateSize - 3).cwiseAbs().maxCoeff() < 1e-15);
}

TEST_CASE("measurements.lander_crater") {
  LanderCraterMeasurement m;
  m.r_crater = Vec3d(300.0, -400.0, 1731200.0);
  m.sigma_rad = 1.0e-3;
  m.config.i_dr = 0;
  m.config.i_dth = 6;

  NavErrorContext nom;
  nom.r = Vec3d(0.0, 0.0, 1730000.0);
  nom.R_b2n = SmallRotation();
  nom.error_state_size = kLanderNavErrorStateSize;

  MatXd H;
  MeasData md = m.Compute(nom, &H);

  SECTION("predicted body line-of-sight is a unit vector consistent with the helper") {
    REQUIRE(md.value.size() == 3);
    Vec3d pred = m.PredictedLosBody(nom.r, nom.R_b2n);
    REQUIRE((md.value - pred).cwiseAbs().maxCoeff() < 1e-12);
    REQUIRE_THAT(md.value.norm(), WithinAbs(1.0, 1e-12));  // unit LOS
  }

  SECTION("noise covariance is sigma_rad^2 * I") {
    REQUIRE(md.covariance.rows() == 3);
    REQUIRE_THAT(md.covariance(0, 0), WithinAbs(1.0e-6, 1e-15));
    REQUIRE_THAT(md.covariance(1, 1), WithinAbs(1.0e-6, 1e-15));
    REQUIRE_THAT(md.covariance(0, 1), WithinAbs(0.0, 1e-18));
  }

  SECTION("position Jacobian block matches finite differences") {
    REQUIRE(H.rows() == 3);
    REQUIRE(H.cols() == kLanderNavErrorStateSize);
    const double eps = 1e-2;
    for (int j = 0; j < 3; ++j) {
      NavErrorContext np = nom, nm = nom;
      np.r(j) += eps;
      nm.r(j) -= eps;
      Vec3d zp = np.R_b2n.transpose() * (m.r_crater - np.r).normalized();
      Vec3d zm = nm.R_b2n.transpose() * (m.r_crater - nm.r).normalized();
      Vec3d fd = (zp - zm) / (2.0 * eps);
      for (int i = 0; i < 3; ++i) REQUIRE_THAT(H(i, j), WithinAbs(fd(i), 1e-6));
    }
  }

  SECTION("degenerate (zero-range) geometry returns an empty value") {
    LanderCraterMeasurement deg = m;
    deg.r_crater = nom.r;  // crater coincident with the lander
    MeasData d = deg.Compute(nom);
    REQUIRE(d.value.size() == 0);
  }
}
