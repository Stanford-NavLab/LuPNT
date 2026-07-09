#include <lupnt/numerics/filters/ekf.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("filters.ekf") {
  Config config = YAML::Load("name: ekf\noutlier_threshold: 100.0\n");
  EKF ekf(config);
  ekf.SetTime(0.0);
  ekf.SetState(State(Vec1(0.0)));
  ekf.SetCovariance(Mat1d::Identity());
  ekf.SetDynamicsFunction([](const State& x, Real t0, Real tf, const State* u, MatXd* F) {
    (void)u;
    if (F) *F = Mat1d::Identity();
    return State(Vec1(x(0) + (tf - t0)));
  });
  ekf.SetProcessNoiseFunction([](const State& x, Real t0, Real tf) {
    (void)x;
    (void)t0;
    (void)tf;
    return Mat1d::Zero();
  });
  ekf.SetMeasurementFunction([](const State& x, MatXd* H, MatXd* R) {
    *H = Mat1d::Identity();
    *R = Mat1d::Identity();
    return Vec1(x(0));
  });

  ekf.Predict(1.0);
  REQUIRE_THAT(ekf.GetState()(0).val(), WithinAbs(1.0, epsilon));
  REQUIRE_THAT(ekf.GetStateJacobian()(0, 0), WithinAbs(1.0, epsilon));

  ekf.Update(Vec1(2.0));
  REQUIRE_THAT(ekf.GetState()(0).val(), WithinAbs(1.5, epsilon));
  REQUIRE_THAT(ekf.GetCovariance()(0, 0), WithinAbs(0.5, epsilon));
  REQUIRE_THAT(ekf.GetMeasurementResidual()(0), WithinAbs(1.0, epsilon));
}
