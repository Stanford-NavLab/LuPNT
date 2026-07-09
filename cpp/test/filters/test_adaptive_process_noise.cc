#include <lupnt/numerics/filters/adaptive_process_noise.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("filters.adaptive_process_noise") {
  ProcessNoise model;
  Mat3d Qa = 0.01 * Mat3d::Identity();
  model.SetAlgorithm(ProcessNoiseAlgorithm::SNC);
  model.SetProcessNoise(Qa);
  model.SetTimeStep(2.0);

  MatXd acc_noise = model.ComputeAccNoise();
  REQUIRE(acc_noise.isApprox(Qa, epsilon));

  MatXd Q = model.ComputeProcessNoise();
  REQUIRE(Q.rows() == 6);
  REQUIRE(Q.cols() == 6);
  REQUIRE_THAT(Q(0, 0), WithinAbs(0.01 * 8.0 / 3.0, epsilon));

  auto proc = model.GetProcessNoiseFunction();
  REQUIRE(proc(State(Vec6::Zero()), 0.0, 2.0).isApprox(Q, epsilon));

  model.SetAlgorithm(ProcessNoiseAlgorithm::ASNC);
  model.SetWaitSteps(0);
  model.SetWindowSize(1);
  model.SetNoiseLimits(0.01, 10.0);
  model.Update(Vec6d::Zero(), Mat6d::Identity(), Mat6d::Identity(), Mat6d::Identity());
  REQUIRE(model.ComputeAccNoise().rows() == 3);
}
