#include <lupnt/filters/ukf.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("filters.ukf") {
  UKF ukf;
  State x0(Vec1(0.0));
  Mat1d P0 = Mat1d::Identity();
  ukf.Initialize(0.0, x0, P0);
  ukf.SetDynamicsFunction([](const State& x, Real t0, Real tf, const State* u, MatXd* F) {
    (void)u;
    if (F) *F = Mat1d::Identity();
    return State(Vec1(x(0) + (tf - t0)));
  });
  ukf.SetProcessNoiseFunction([](const State& x, Real t0, Real tf) {
    (void)x;
    (void)t0;
    (void)tf;
    return Mat1d::Zero();
  });
  ukf.SetMeasurementFunction([](const State& x, MatXd* H, MatXd* R) {
    *H = Mat1d::Identity();
    *R = Mat1d::Identity();
    return Vec1(x(0));
  });

  MatX sigma = ukf.ComputeSigmaPoints(x0, P0);
  REQUIRE(sigma.rows() == 1);
  REQUIRE(sigma.cols() == 3);
  REQUIRE_THAT(sigma(0, 0).val(), WithinAbs(0.0, epsilon));

  ukf.Predict(1.0);
  REQUIRE_THAT(ukf.GetState()(0).val(), WithinAbs(1.0, 1e-9));
  ukf.Update(Vec1(2.0));
  REQUIRE_THAT(ukf.GetState()(0).val(), WithinAbs(1.5, 1e-6));
  REQUIRE_THAT(ukf.GetCovariance()(0, 0), WithinAbs(0.5, 1e-6));
}
