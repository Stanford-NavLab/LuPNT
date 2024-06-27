#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch::Matchers;
using namespace lupnt;

static void RequireNearRealVec(const VecX& a, const VecX& b, double abs_error) {
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.size(); ++i) {
    REQUIRE_THAT(a[i].val(), WithinAbs(b[i].val(), abs_error));
  }
}

static void RequireNearRealMat(const MatX& a, const MatX& b, double abs_error) {
  REQUIRE(a.rows() == b.rows());
  REQUIRE(a.cols() == b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      REQUIRE_THAT(a(i, j).val(), WithinAbs(b(i, j).val(), abs_error));
    }
  }
}

static void RequireNearDoubleVec(const VecXd& a, const VecXd& b,
                                 double abs_error) {
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.size(); ++i) {
    REQUIRE_THAT(a[i], WithinAbs(b[i], abs_error));
  }
}

static void RequireNearDoubleMat(const VecXd& a, const VecXd& b,
                                 double abs_error) {
  REQUIRE(a.rows() == b.rows());
  REQUIRE(a.cols() == b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      REQUIRE_THAT(a(i, j), WithinAbs(b(i, j), abs_error));
    }
  }
}

static void RequireNearReal(Real a, Real b, double abs_error) {
  REQUIRE_THAT(a.val(), WithinAbs(b.val(), abs_error));
}

static void RequireNearDouble(double a, double b, double abs_error) {
  REQUIRE_THAT(a, WithinAbs(b, abs_error));
}

// inputs: vector, function(vector), jacobian
static void NumericalJacobian(
    std::function<void(VecX&, Real)> propagate_function, const VecX& vec,
    Real dt, Mat6d& jacobian, double eps = 1e-6) {
  int n = vec.size();
  VecX vec_p;
  VecX vec_m;

  for (int i = 0; i < n; i++) {
    vec_p = vec;
    vec_m = vec;
    vec_p(i) += eps;
    vec_m(i) -= eps;

    propagate_function(vec_p, dt);
    propagate_function(vec_m, dt);

    jacobian.col(i) = (vec_p - vec_m).cast<double>() / (2.0 * eps);
  }
}