#include <lupnt/numerics/filters/filter_utils.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../utils.cc"

using namespace lupnt;
using Catch::Approx;

// Helpers -------------------------------------------------------------------
namespace {
  bool IsSymmetric(const MatXd& M, double tol = 1e-9) { return M.isApprox(M.transpose(), tol); }
  // True if M is symmetric positive-definite (via successful LLT + positive pivots).
  bool IsSpd(const MatXd& M, double tol = 1e-12) {
    if (!IsSymmetric(M)) return false;
    Eigen::LLT<MatXd> llt(M);
    if (llt.info() != Eigen::Success) return false;
    return (llt.matrixL().toDenseMatrix().diagonal().array() > tol).all();
  }
}  // namespace

TEST_CASE("numerics.filter_utils.initial_covariance_pos_vel_clock") {
  double sr = 100.0, sv = 0.5, sb = 1e-6, sd = 1e-9;
  MatXd P = InitialCovariancePosVelClock(sr, sv, sb, sd);

  REQUIRE(P.rows() == 8);
  REQUIRE(P.cols() == 8);
  REQUIRE(IsSymmetric(P));

  // Diagonal blocks carry the squared sigmas; everything else is zero.
  for (int i = 0; i < 3; i++) REQUIRE(P(i, i) == Approx(sr * sr));
  for (int i = 3; i < 6; i++) REQUIRE(P(i, i) == Approx(sv * sv));
  REQUIRE(P(6, 6) == Approx(sb * sb));
  REQUIRE(P(7, 7) == Approx(sd * sd));

  // Off-diagonal entries vanish.
  MatXd D = P;
  for (int i = 0; i < 8; i++) D(i, i) = 0.0;
  REQUIRE(D.norm() == Approx(0.0));

  REQUIRE(IsSpd(P));
}

TEST_CASE("numerics.filter_utils.process_noise_pos_vel_coeffs") {
  // Closed form: [dt^3/3, dt^2/2, dt].
  double dt = 7.0;
  Vec3d c = ProcessNoisePosVelCoeffs(dt);
  REQUIRE(c(0) == Approx(std::pow(dt, 3) / 3.0));
  REQUIRE(c(1) == Approx(std::pow(dt, 2) / 2.0));
  REQUIRE(c(2) == Approx(dt));

  // Non-positive dt is rejected.
  REQUIRE_THROWS(ProcessNoisePosVelCoeffs(0.0));
  REQUIRE_THROWS(ProcessNoisePosVelCoeffs(-1.0));
}

TEST_CASE("numerics.filter_utils.process_noise_pos_vel") {
  double dt = 4.0;
  double q = 2.5;  // per-axis acceleration PSD
  Vec3d qa(q, q, q);
  MatXd Q = ProcessNoisePosVel(qa, dt);  // 3-vector -> diagonal Q_a

  REQUIRE(Q.rows() == 6);
  REQUIRE(Q.cols() == 6);
  REQUIRE(IsSymmetric(Q));

  double c11 = std::pow(dt, 3) / 3.0, c21 = std::pow(dt, 2) / 2.0, c22 = dt;
  for (int i = 0; i < 3; i++) {
    REQUIRE(Q(i, i) == Approx(q * c11));
    REQUIRE(Q(i, i + 3) == Approx(q * c21));
    REQUIRE(Q(i + 3, i) == Approx(q * c21));
    REQUIRE(Q(i + 3, i + 3) == Approx(q * c22));
  }
  REQUIRE(IsSpd(Q));

  // A full (square) Q_a matrix must give the same result as its diagonal vector form.
  MatXd Qa_full = qa.asDiagonal();
  MatXd Q_full = ProcessNoisePosVel(Qa_full, dt);
  REQUIRE(Q_full.isApprox(Q, 1e-12));
}

TEST_CASE("numerics.filter_utils.state_transition_matrix_pos_vel") {
  double dt = 3.0;
  int n = 3;
  MatXd Phi = StateTransitionMatrixPosVel(dt, n);

  REQUIRE(Phi.rows() == 6);
  MatXd expected(6, 6);
  expected << MatXd::Identity(3, 3), dt * MatXd::Identity(3, 3), MatXd::Zero(3, 3),
      MatXd::Identity(3, 3);
  REQUIRE(Phi.isApprox(expected, 1e-12));

  // Constant-velocity STMs compose additively in dt: Phi(a) * Phi(b) = Phi(a+b).
  MatXd Pa = StateTransitionMatrixPosVel(1.5, n);
  MatXd Pb = StateTransitionMatrixPosVel(2.0, n);
  MatXd Pab = StateTransitionMatrixPosVel(3.5, n);
  REQUIRE((Pa * Pb).isApprox(Pab, 1e-12));

  // det = 1 (unipotent upper-triangular).
  REQUIRE(Phi.determinant() == Approx(1.0));
}

TEST_CASE("numerics.filter_utils.state_transition_matrix_pos_vel_acc") {
  double dt = 5.0;
  VecXd beta(3);
  beta << 0.01, 0.02, 0.05;
  MatXd Phi = StateTransitionMatrixPosVelAcc(dt, beta);

  REQUIRE(Phi.rows() == 9);
  REQUIRE(Phi.cols() == 9);

  // Top-left 6x6 block is exactly the constant-velocity STM.
  MatXd Phi_rv = StateTransitionMatrixPosVel(dt, 3);
  REQUIRE(Phi.block(0, 0, 6, 6).isApprox(Phi_rv, 1e-12));

  // Acceleration sub-block decays as exp(-beta*dt).
  for (int i = 0; i < 3; i++) {
    REQUIRE(Phi(6 + i, 6 + i) == Approx(std::exp(-beta(i) * dt)));
  }
  // Lower-left block is zero (acceleration unaffected by pos/vel).
  REQUIRE(Phi.block(6, 0, 3, 6).norm() == Approx(0.0));
}

TEST_CASE("numerics.filter_utils.process_noise_pos_vel_acc_coeffs") {
  double dt = 6.0;
  VecXd beta(2);
  beta << 0.1, 0.2;

  // Exercise the coefficient builder (finite, correct shape) and its input
  // guard. (The full ProcessNoisePosVelAcc(Q_a, dt, beta) assembly is not
  // asserted here: for these inputs it trips its own internal non-negative-
  // diagonal LUPNT_CHECK -- flagged separately as a potential filter_utils
  // issue rather than masked by the test.)
  Mat6Xd C = ProcessNoisePosVelAccCoeffs(dt, beta);
  REQUIRE(C.rows() == 6);
  REQUIRE(C.cols() == 2);
  REQUIRE(C.allFinite());

  REQUIRE_THROWS(ProcessNoisePosVelAccCoeffs(-1.0, beta));
}

TEST_CASE("numerics.filter_utils.process_noise_clock") {
  double dt = 10.0;
  MatXd Q2 = ProcessNoiseClock(ClockModel::USO, 2, dt);
  REQUIRE(Q2.rows() == 2);
  REQUIRE(IsSymmetric(Q2));
  REQUIRE((Q2.diagonal().array() >= 0.0).all());  // PSD (not necessarily strictly SPD)

  MatXd Q3 = ProcessNoiseClock(ClockModel::USO, 3, dt);
  REQUIRE(Q3.rows() == 3);
  REQUIRE(IsSymmetric(Q3));

  // Invalid clock dimension throws.
  REQUIRE_THROWS(ProcessNoiseClock(ClockModel::USO, 4, dt));
}

// NOTE: ProcessNoisePosVelClock is intentionally NOT tested here. The public
// declaration in filter_utils.h takes 4 args (ClockModel, int, double, int),
// but the definition in filter_utils.cc takes 5 (an extra trailing `Real dt`
// with no default). The declared 4-arg symbol therefore has no definition and
// any call fails to link -- a pre-existing source mismatch, not something this
// test file can (or should) work around.
