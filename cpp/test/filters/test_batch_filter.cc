#include <lupnt/filters/batch_filter.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("filters.batch_filter") {
  MatXd H(2, 1);
  H << 1.0, 1.0;
  VecXd residuals(2);
  residuals << 2.0, 4.0;
  VecXd weights = VecXd::Ones(2);
  VecXd dx = SolveWeightedLeastSquares(H, residuals, weights);
  REQUIRE_THAT(dx(0), WithinAbs(3.0, epsilon));

  VecXd initial(1);
  initial << 0.0;
  VecXd current(1);
  current << 1.0;
  auto [H_aug, r_aug, w_aug]
      = AddInitializationConstraints(H, residuals, weights, initial, current, VecXd::Ones(1));
  REQUIRE(H_aug.rows() == 3);
  REQUIRE_THAT(r_aug(2), WithinAbs(-1.0, epsilon));
  REQUIRE_THAT(w_aug(2), WithinAbs(1.0, epsilon));

  BatchFilterConfig config;
  config.use_initialization = false;
  config.convergence_tol = 1e-12;
  std::vector<VecXd> measurements{Vec1(2.0), Vec1(4.0)};
  std::vector<VecXd> measurement_weights{Vec1(1.0), Vec1(1.0)};
  auto model = [](const VecXd& state, int meas_idx) {
    (void)meas_idx;
    MatXd H_i = MatXd::Ones(1, 1);
    return std::make_pair(Vec1(state(0)), H_i);
  };
  BatchFilterResults results = RunBatchFilter(initial, Vec1(1.0), measurements, measurement_weights,
                                              model, config, Vec1(3.0));
  REQUIRE(results.converged);
  REQUIRE_THAT(results.state_estimate(0), WithinAbs(3.0, epsilon));
  REQUIRE_THAT(results.state_errors(0), WithinAbs(0.0, epsilon));
  REQUIRE(std::isinf(CalculatePDOP(Mat1d::Identity())));
}
