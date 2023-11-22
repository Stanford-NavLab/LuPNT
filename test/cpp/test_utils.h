#pragma once

#include <gtest/gtest.h>
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

using namespace lupnt;

static void EXPECT_NEAR_ADVEC(const VectorXreal& a, const VectorXreal& b,
                              double abs_error) {
  EXPECT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); ++i) {
    EXPECT_NEAR(a[i].val(), b[i].val(), abs_error);
  }
}

static void EXPECT_NEAR_ADMAT(const MatrixXreal& a, const MatrixXreal& b,
                              double abs_error) {
  EXPECT_EQ(a.rows(), b.rows());
  EXPECT_EQ(a.cols(), b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      EXPECT_NEAR(a(i, j).val(), b(i, j).val(), abs_error);
    }
  }
}

static void EXPECT_NEAR_EIGENVEC(const VectorXd& a, const VectorXd& b,
                                 double abs_error) {
  EXPECT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); ++i) {
    EXPECT_NEAR(a[i], b[i], abs_error);
  }
}

static void EXPECT_NEAR_EIGENMAT(const MatrixXd& a, const MatrixXd& b,
                                 double abs_error) {
  EXPECT_EQ(a.rows(), b.rows());
  EXPECT_EQ(a.cols(), b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      EXPECT_NEAR(a(i, j), b(i, j), abs_error);
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