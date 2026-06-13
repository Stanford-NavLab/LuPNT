#include <lupnt/filters/filter.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

namespace {
  class TestFilter : public Filter {
  public:
    using Filter::Filter;
    void Predict(Real t, const State* u = nullptr) override {
      (void)u;
      t_ = t;
      x_prior_ = x_;
      P_prior_ = P_;
    }
    void Update(const VecX& z_true) override {
      (void)z_true;
      x_post_ = x_;
      P_post_ = P_;
    }
  };
}  // namespace

TEST_CASE("filters.filter") {
  Config config = YAML::Load("name: base_filter\n");
  TestFilter filter(config);
  State x(Vec2(1.0, 2.0));
  Mat2d P = 4.0 * Mat2d::Identity();

  filter.SetTime(5.0);
  filter.SetState(x);
  filter.SetCovariance(P);
  filter.SetMeasurementFunction([](const State& state, MatXd* H, MatXd* R) {
    *H = MatXd::Zero(1, state.size());
    (*H)(0, 0) = 1.0;
    *R = MatXd::Identity(1, 1) * 4.0;
    return Vec1(state(0));
  });

  REQUIRE(filter.GetName() == "base_filter");
  REQUIRE_THAT(filter.GetTime().val(), WithinAbs(5.0, epsilon));
  REQUIRE(filter.GetState().isApprox(x, epsilon));
  REQUIRE(filter.GetCovariance().isApprox(P, epsilon));
  REQUIRE_THAT(filter.ComputeResidualRMS(Vec1(3.0), x), WithinAbs(1.0, epsilon));
  REQUIRE_THAT(filter.ComputeResidualRMS(VecXd(), x), WithinAbs(0.0, epsilon));
}
