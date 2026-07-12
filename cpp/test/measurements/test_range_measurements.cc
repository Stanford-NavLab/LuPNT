// Unit tests for the two geometry-only measurement models that carry closed-form
// / autodiff design matrices: `GroundStationRangeMeasurement` (used by the ex7
// ground-station batch OD) and `IslCrosslinkMeasurement` (used by the ex8
// distributed ISL ODTS). Both are self-contained: a hand-computed state gives a
// known observable vector, covariance, and Jacobian, so no data files are needed.

#include <lupnt/measurements/crosslink_measurement.h>
#include <lupnt/measurements/ground_range_measurement.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Numerical Jacobian of a measurement model's value w.r.t. the state vector,
  // by central differences. Returns an [n_z x n_x] matrix.
  MatXd NumericalMeasJacobian(const Measurement& meas, const State& x0, double eps = 1e-4) {
    const int nx = static_cast<int>(x0.size());
    const int nz = static_cast<int>(meas.Compute(x0).value.size());
    MatXd J(nz, nx);
    for (int j = 0; j < nx; ++j) {
      State xp = x0, xm = x0;
      xp(j) += eps;
      xm(j) -= eps;
      VecXd zp = meas.Compute(xp).value;
      VecXd zm = meas.Compute(xm).value;
      J.col(j) = (zp - zm) / (2.0 * eps);
    }
    return J;
  }
}  // namespace

TEST_CASE("measurements.ground_range.value_and_jacobian") {
  // Station at the origin (at rest); target on the 3-4-5 triangle with a velocity
  // aligned so the range rate is exactly the closing speed.
  GroundStationRangeMeasurement::Config cfg;
  cfg.use_range = true;
  cfg.use_range_rate = true;
  cfg.reference_state = Vec6d::Zero();
  cfg.range_sigma_m = 10.0;
  cfg.range_rate_sigma_mps = 1.0e-3;
  GroundStationRangeMeasurement meas(cfg);

  Cart6 x(Vec6(3.0, 4.0, 0.0, 0.6, 0.8, 0.0), Frame::MOON_CI);

  MatXd H;
  MeasData md = meas.Compute(x, &H);

  SECTION("observable vector") {
    REQUIRE(md.value.size() == 2);
    REQUIRE_THAT(md.value(0), WithinAbs(5.0, 1e-9));  // range = |dr| = 5
    REQUIRE_THAT(md.value(1), WithinAbs(1.0, 1e-9));  // range rate = dr.dv/|dr| = 1
  }

  SECTION("noise covariance is the configured diagonal") {
    REQUIRE(md.covariance.rows() == 2);
    REQUIRE_THAT(md.covariance(0, 0), WithinAbs(100.0, 1e-12));
    REQUIRE_THAT(md.covariance(1, 1), WithinAbs(1.0e-6, 1e-18));
    REQUIRE_THAT(md.covariance(0, 1), WithinAbs(0.0, 1e-18));
  }

  SECTION("closed-form design matrix") {
    REQUIRE(H.rows() == 2);
    REQUIRE(H.cols() == 6);
    // range row = [u^T, 0], u = dr/|dr| = (0.6, 0.8, 0)
    REQUIRE_THAT(H(0, 0), WithinAbs(0.6, 1e-9));
    REQUIRE_THAT(H(0, 1), WithinAbs(0.8, 1e-9));
    REQUIRE_THAT(H(0, 3), WithinAbs(0.0, 1e-9));
    // range-rate velocity partial = u; position partial = (dv - rdot*u)/rho = 0 here
    REQUIRE_THAT(H(1, 0), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(H(1, 3), WithinAbs(0.6, 1e-9));
    REQUIRE_THAT(H(1, 4), WithinAbs(0.8, 1e-9));
  }

  SECTION("closed-form matches numerical differentiation") {
    MatXd Hnum = NumericalMeasJacobian(meas, x);
    REQUIRE((H - Hnum).cwiseAbs().maxCoeff() < 1e-5);
  }

  SECTION("range-only config drops the range-rate row") {
    GroundStationRangeMeasurement::Config c = cfg;
    c.use_range_rate = false;
    GroundStationRangeMeasurement m(c);
    MeasData d = m.Compute(x);
    REQUIRE(d.value.size() == 1);
    REQUIRE_THAT(d.value(0), WithinAbs(5.0, 1e-9));
  }
}

TEST_CASE("measurements.isl_crosslink.value_covariance_jacobian") {
  // Hub block (index 0) at rest at the origin; one linked satellite on the 3-4-5
  // triangle with a closing velocity. Each block is [r(3), v(3), cb, cd] (size 8).
  IslCrosslinkMeasurement::Config cfg;
  cfg.sub_state_size = 8;
  cfg.n_links = 1;
  cfg.sigma_range_m = 2.0;
  cfg.sigma_range_rate_mps = 1.0e-3;
  cfg.sigma_pseudorange_m = 50.0;
  cfg.anchor_pos_mci = {Vec3d(0.0, 0.0, 10.0)};  // one known-position anchor

  IslCrosslinkMeasurement meas(cfg);

  State x(16);
  x.setZero();
  // hub: origin, at rest, clock bias 1e-6 s
  x(6) = 1.0e-6;
  // linked satellite block (offset 8): r = (3,4,0), v = (0.6, 0.8, 0)
  x.segment(8, 3) = Vec3(3.0, 4.0, 0.0);
  x.segment(11, 3) = Vec3(0.6, 0.8, 0.0);

  MatXd H;
  MeasData md = meas.Compute(x, &H);

  SECTION("row layout: 2*n_links crosslink rows + one anchor pseudorange") {
    REQUIRE(md.value.size() == 3);
    REQUIRE_THAT(md.value(0), WithinAbs(5.0, 1e-6));  // crosslink range
    REQUIRE_THAT(md.value(1), WithinAbs(1.0, 1e-6));  // crosslink range rate
    // anchor pseudorange = |r_hub - r_anchor| + C*b_hub = 10 + C*1e-6
    REQUIRE_THAT(md.value(2), WithinAbs(10.0 + C * 1.0e-6, 1e-3));
  }

  SECTION("diagonal noise covariance from the configured sigmas") {
    REQUIRE(md.covariance.rows() == 3);
    REQUIRE_THAT(md.covariance(0, 0), WithinAbs(4.0, 1e-9));
    REQUIRE_THAT(md.covariance(1, 1), WithinAbs(1.0e-6, 1e-15));
    REQUIRE_THAT(md.covariance(2, 2), WithinAbs(2500.0, 1e-6));
  }

  SECTION("autodiff Jacobian matches numerical differentiation") {
    REQUIRE(H.rows() == 3);
    REQUIRE(H.cols() == 16);
    MatXd Hnum = NumericalMeasJacobian(meas, x, 1e-3);
    REQUIRE((H - Hnum).cwiseAbs().maxCoeff() < 1e-3);
  }

  SECTION("optional time/frequency-transfer rows are appended") {
    IslCrosslinkMeasurement::Config c = cfg;
    c.include_time_transfer = true;
    c.include_frequency_transfer = true;
    IslCrosslinkMeasurement m(c);
    State x2 = x;
    x2(8 + 6) = 2.0e-6;  // linked clock bias
    x2(8 + 7) = 3.0e-9;  // linked clock drift
    MeasData d = m.Compute(x2);
    // 2 crosslink + 1 anchor + 1 time-transfer + 1 frequency-transfer
    REQUIRE(d.value.size() == 5);
    // time transfer = C*(b_hub - b_link) = C*(1e-6 - 2e-6)
    REQUIRE_THAT(d.value(3), WithinAbs(C * (1.0e-6 - 2.0e-6), 1e-6));
    // frequency transfer = C*(d_hub - d_link) = C*(0 - 3e-9)
    REQUIRE_THAT(d.value(4), WithinAbs(C * (0.0 - 3.0e-9), 1e-9));
  }
}
