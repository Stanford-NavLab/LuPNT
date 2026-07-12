#include <lupnt/lupnt.h>
#include <lupnt/numerics/filters/ekf.h>
#include <lupnt/numerics/filters/udu_filter.h>
#include <lupnt/numerics/filters/udu_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  // A generic symmetric-positive-definite matrix built as A^T A + n I so it is well
  // conditioned and has non-trivial off-diagonal structure.
  MatXd MakeSpd(int n, double seed) {
    MatXd A(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) A(i, j) = std::sin(seed + i * 1.7 + j * 0.9) + 0.3 * (i - j);
    return A.transpose() * A + n * MatXd::Identity(n, n);
  }

  // Constant-velocity dynamics [pos; vel] with STM [[1, dt], [0, 1]].
  FilterDynamicsFunction ConstVelDynamics() {
    return [](const State& x, Real t0, Real tf, const State*, MatXd* F) {
      double dt = (tf - t0).val();
      MatXd Phi(2, 2);
      Phi << 1.0, dt, 0.0, 1.0;
      if (F != nullptr) *F = Phi;
      State x_next(2);
      x_next << x(0) + dt * x(1), x(1);
      return x_next;
    };
  }

  // Diagonal (position/velocity) process noise, as required by the UDU time update.
  ProcessNoiseFunction DiagProcessNoise(double q_pos, double q_vel) {
    return [q_pos, q_vel](const State&, Real, Real) {
      MatXd Q = MatXd::Zero(2, 2);
      Q(0, 0) = q_pos;
      Q(1, 1) = q_vel;
      return Q;
    };
  }

  FilterMeasurementFunction PositionMeasurement(double r) {
    return [r](const State& x, MatXd* H_out, MatXd* R_out) {
      if (H_out != nullptr) {
        *H_out = MatXd::Zero(1, x.size());
        (*H_out)(0, 0) = 1.0;
      }
      if (R_out != nullptr) *R_out = r * MatXd::Identity(1, 1);
      VecXd y(1);
      y(0) = x(0).val();
      return y;
    };
  }

}  // namespace

TEST_CASE("numerics.udu_utils.decomposition_roundtrip") {
  SECTION("general SPD matrix reconstructs exactly") {
    for (int n : {1, 2, 4, 6}) {
      MatXd P = MakeSpd(n, 0.5 * n);
      auto [D, U] = UDUDecomposition(P);
      // U is unit upper triangular
      REQUIRE(U.rows() == n);
      REQUIRE(D.size() == n);
      for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(U(i, i), WithinAbs(1.0, 1.0e-12));
        for (int j = 0; j < i; ++j) REQUIRE_THAT(U(i, j), WithinAbs(0.0, 1.0e-12));
        REQUIRE(D(i) > 0.0);  // SPD -> strictly positive diagonal factor
      }
      MatXd P_rec = UDUReconstruct(U, D);
      REQUIRE_THAT((P - P_rec).norm(), WithinAbs(0.0, 1.0e-9));
    }
  }

  SECTION("diagonal covariance factors to identity U and D = diag(P)") {
    MatXd P = MatXd::Zero(3, 3);
    P.diagonal() << 2.0, 5.0, 0.25;
    auto [D, U] = UDUDecomposition(P);
    REQUIRE_THAT((U - MatXd::Identity(3, 3)).norm(), WithinAbs(0.0, 1.0e-12));
    REQUIRE_THAT(D(0), WithinAbs(2.0, 1.0e-12));
    REQUIRE_THAT(D(1), WithinAbs(5.0, 1.0e-12));
    REQUIRE_THAT(D(2), WithinAbs(0.25, 1.0e-12));
  }

  SECTION("numerically-zero variance is clamped to exactly zero") {
    MatXd P = MatXd::Zero(2, 2);
    P(0, 0) = 3.0;
    P(1, 1) = 1.0e-20;  // below the 1e-16 clamp threshold
    auto [D, U] = UDUDecomposition(P);
    REQUIRE(D(1) == 0.0);
  }
}

TEST_CASE("numerics.udu_utils.modified_gram_schmidt") {
  // MWGS must satisfy U diag(D) U^T = Y D_tilde Y^T without ever forming the product.
  const int n = 3;
  const int m = 2;
  MatXd Y(n, n + m);
  Y << 1.0, 0.2, -0.3, 0.5, 0.1, 0.0, 1.0, 0.4, -0.2, 0.3, 0.1, -0.1, 1.0, 0.2, 0.6;
  VecXd d_tilde_diag(n + m);
  d_tilde_diag << 3.0, 1.5, 0.8, 0.5, 0.2;
  MatXd D_tilde = d_tilde_diag.asDiagonal();

  auto [D, U] = ModifiedGramSchmidt(D_tilde, Y);
  MatXd lhs = UDUReconstruct(U, D);
  MatXd rhs = Y * D_tilde * Y.transpose();
  REQUIRE_THAT((lhs - rhs).norm(), WithinAbs(0.0, 1.0e-10));

  // U unit upper triangular, D non-negative
  for (int i = 0; i < n; ++i) {
    REQUIRE_THAT(U(i, i), WithinAbs(1.0, 1.0e-12));
    REQUIRE(D(i) >= 0.0);
  }
}

TEST_CASE("numerics.udu_utils.solve") {
  MatXd P = MakeSpd(4, 2.0);
  auto [D, U] = UDUDecomposition(P);

  SECTION("solves P X = B for a matrix RHS") {
    MatXd B(4, 2);
    B << 1.0, -2.0, 0.5, 3.0, -1.5, 0.25, 2.0, -0.75;
    MatXd X = UDUSolve(U, D, B);
    REQUIRE_THAT((P * X - B).norm(), WithinAbs(0.0, 1.0e-8));
  }

  SECTION("recovers the identity when B = P") {
    MatXd X = UDUSolve(U, D, P);
    REQUIRE_THAT((X - MatXd::Identity(4, 4)).norm(), WithinAbs(0.0, 1.0e-8));
  }
}

TEST_CASE("numerics.udu_filter.matches_ekf_over_multiple_cycles") {
  // A constant-velocity target tracked by both a plain EKF and the UDU EKF, with identical
  // (diagonal) process noise and a scalar position measurement. Their posterior mean and
  // covariance must agree at every epoch to tight tolerance.
  State x0(2);
  x0 << 0.0, 1.0;
  MatXd P0(2, 2);
  P0 << 10.0, 0.0, 0.0, 4.0;
  const double q_pos = 0.01, q_vel = 0.02, r = 0.5;

  EKF ekf;
  ekf.SetState(x0);
  ekf.SetCovariance(P0);
  ekf.SetDynamicsFunction(ConstVelDynamics());
  ekf.SetProcessNoiseFunction(DiagProcessNoise(q_pos, q_vel));
  ekf.SetMeasurementFunction(PositionMeasurement(r));

  UDUEKF udu;
  udu.SetState(x0);
  udu.SetCovariance(P0);
  udu.SetDynamicsFunction(ConstVelDynamics());
  udu.SetProcessNoiseFunction(DiagProcessNoise(q_pos, q_vel));
  udu.SetMeasurementFunction(PositionMeasurement(r));

  const std::vector<double> measurements = {0.9, 2.1, 2.8, 4.2, 4.9};
  for (int k = 0; k < static_cast<int>(measurements.size()); ++k) {
    Real t = static_cast<double>(k + 1);
    ekf.Predict(t);
    udu.Predict(t);

    // Predicted covariance also agrees.
    REQUIRE_THAT((ekf.GetCovariance() - udu.GetCovariance()).norm(), WithinAbs(0.0, 1.0e-9));

    VecXd z(1);
    z(0) = measurements[k];
    ekf.Update(z);
    udu.Update(z);

    REQUIRE_THAT((ekf.GetState() - udu.GetState()).norm().val(), WithinAbs(0.0, 1.0e-9));
    REQUIRE_THAT((ekf.GetCovariance() - udu.GetCovariance()).norm(), WithinAbs(0.0, 1.0e-9));

    // The UDU covariance must stay symmetric positive definite.
    MatXd P = udu.GetCovariance();
    REQUIRE_THAT((P - P.transpose()).norm(), WithinAbs(0.0, 1.0e-10));
    Eigen::SelfAdjointEigenSolver<MatXd> es(P);
    REQUIRE(es.eigenvalues().minCoeff() > 0.0);

    // The UDU factors reconstruct the reported covariance.
    REQUIRE_THAT((UDUReconstruct(udu.GetUFactor(), udu.GetDFactor()) - P).norm(),
                 WithinAbs(0.0, 1.0e-10));
  }
}

TEST_CASE("numerics.udu_filter.update_reduces_position_variance") {
  State x0(2);
  x0 << 0.0, 0.0;
  MatXd P0(2, 2);
  P0 << 9.0, 0.0, 0.0, 4.0;

  UDUEKF udu;
  udu.SetState(x0);
  udu.SetCovariance(P0);
  udu.SetDynamicsFunction(ConstVelDynamics());
  udu.SetProcessNoiseFunction(DiagProcessNoise(0.0, 0.0));
  udu.SetMeasurementFunction(PositionMeasurement(1.0));

  udu.Predict(1.0);
  double var_pos_prior = udu.GetCovariance()(0, 0);
  VecXd z(1);
  z(0) = 0.5;
  udu.Update(z);
  double var_pos_post = udu.GetCovariance()(0, 0);

  REQUIRE(var_pos_post < var_pos_prior);  // measuring position shrinks its variance
  REQUIRE(var_pos_post > 0.0);
}

TEST_CASE("numerics.udu_filter.stochastic_cloning_predict") {
  // The cloned filter carries [x_current; x_previous]. After a predict the previous block
  // equals the prior current state and the current block is the propagated state.
  const int base_n = 2;
  State x0(2 * base_n);
  x0 << 0.0, 1.0, 0.0, 1.0;  // current = previous = [pos=0, vel=1]
  MatXd P0 = MatXd::Zero(2 * base_n, 2 * base_n);
  P0.diagonal() << 10.0, 4.0, 10.0, 4.0;

  UDUStochasticCloningEKF udu;
  udu.SetBaseStateSize(base_n);
  udu.SetState(x0);
  udu.SetCovariance(P0);
  udu.SetDynamicsFunction(ConstVelDynamics());
  udu.SetProcessNoiseFunction(DiagProcessNoise(0.01, 0.02));

  REQUIRE(udu.GetBaseStateSize() == base_n);

  udu.Predict(1.0);
  State x = udu.GetState();
  REQUIRE(x.size() == 2 * base_n);
  // current block: pos advanced by dt*vel = 1, vel unchanged
  REQUIRE_THAT(x(0).val(), WithinAbs(1.0, 1.0e-9));
  REQUIRE_THAT(x(1).val(), WithinAbs(1.0, 1.0e-9));
  // previous block: the pre-predict current state [0, 1]
  REQUIRE_THAT(x(2).val(), WithinAbs(0.0, 1.0e-9));
  REQUIRE_THAT(x(3).val(), WithinAbs(1.0, 1.0e-9));

  State x_curr = udu.GetCurrentState();
  REQUIRE(x_curr.size() == base_n);
  REQUIRE_THAT(x_curr(0).val(), WithinAbs(1.0, 1.0e-9));

  MatXd P_curr = udu.GetCurrentCovariance();
  REQUIRE(P_curr.rows() == base_n);
  REQUIRE(P_curr.cols() == base_n);
  Eigen::SelfAdjointEigenSolver<MatXd> es(P_curr);
  REQUIRE(es.eigenvalues().minCoeff() > 0.0);
}
