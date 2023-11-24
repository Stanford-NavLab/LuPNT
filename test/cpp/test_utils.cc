#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch::Matchers;
using namespace lupnt;

static void EXPECT_NEAR_ADVEC(const VectorXreal& a, const VectorXreal& b,
                              double abs_error) {
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.size(); ++i) {
    REQUIRE_THAT(a[i].val(), WithinAbs(b[i].val(), abs_error));
  }
}

static void EXPECT_NEAR_ADMAT(const MatrixXreal& a, const MatrixXreal& b,
                              double abs_error) {
  REQUIRE(a.rows() == b.rows());
  REQUIRE(a.cols() == b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      REQUIRE_THAT(a(i, j).val(), WithinAbs(b(i, j).val(), abs_error));
    }
  }
}

static void EXPECT_NEAR_EIGENVEC(const VectorXd& a, const VectorXd& b,
                                 double abs_error) {
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.size(); ++i) {
    REQUIRE_THAT(a[i], WithinAbs(b[i], abs_error));
  }
}

static void EXPECT_NEAR_EIGENMAT(const MatrixXd& a, const MatrixXd& b,
                                 double abs_error) {
  REQUIRE(a.rows() == b.rows());
  REQUIRE(a.cols() == b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      REQUIRE_THAT(a(i, j), WithinAbs(b(i, j), abs_error));
    }
  }
}

// inputs: vector, function(vector), jacobian
static void NumericalJacobian(
    std::function<void(VectorXreal&, real)> propagate_function,
    const VectorXreal& vec, real dt, Matrix6d& jacobian, double eps = 1e-6) {
  int n = vec.size();
  VectorXreal vec_p;
  VectorXreal vec_m;

  for (int i = 0; i < n; i++) {
    vec_p = vec;
    vec_m = vec;
    vec_p(i) += eps;
    vec_m(i) -= eps;

    propagate_function(vec_p, dt);
    propagate_function(vec_m, dt);

    jacobian.col(i) = lupnt::toEigen(vec_p - vec_m) / (2.0 * eps);
  }
}