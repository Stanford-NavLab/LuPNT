#include <lupnt/lupnt.h>
#include <lupnt/numerics/filters/srif.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Linear dynamics x_{k+1} = F x_k with constant state-transition matrix F.
  FilterDynamicsFunction LinearDynamics(const MatXd& F) {
    return [F](const State& x, Real, Real, const State*, MatXd* F_out) -> State {
      if (F_out != nullptr) *F_out = F;
      return State((F * x.cast<double>()).cast<Real>());
    };
  }

  ProcessNoiseFunction ConstProcessNoise(const MatXd& Q) {
    return [Q](const State&, Real, Real) { return Q; };
  }

  FilterMeasurementFunction LinearMeasurement(const MatXd& H, const MatXd& R) {
    return [H, R](const State& x, MatXd* H_out, MatXd* R_out) {
      if (H_out != nullptr) *H_out = H;
      if (R_out != nullptr) *R_out = R;
      return VecXd(H * x.cast<double>());
    };
  }

  // Drive a filter through a forward pass and a backward RTS smoother, collecting the
  // smoothed states/covariances. Works for any KalmanFilter-derived type (EKF, SRIF).
  template <typename F>
  void RunAndSmooth(F& filt, const State& x0, const MatXd& P0, const FilterDynamicsFunction& fd,
                    const ProcessNoiseFunction& fp, const FilterMeasurementFunction& fm,
                    const std::vector<double>& times, const std::vector<VecXd>& zs,
                    std::vector<VecXd>& x_sm, std::vector<MatXd>& P_sm) {
    filt.SetState(x0);
    filt.SetCovariance(P0);
    filt.SetDynamicsFunction(fd);
    filt.SetProcessNoiseFunction(fp);
    filt.SetMeasurementFunction(fm);
    filt.SetOutlierThreshold(1.0e12);  // disable outlier rejection for the comparison
    filt.SetTime(times[0]);
    const int N = static_cast<int>(times.size());
    filt.InitializeLogger(N);
    for (int i = 0; i < N; ++i) {
      if (i > 0) filt.Predict(times[i]);
      filt.Update(zs[i]);
      filt.LogFilterEstimate(i);
    }
    filt.InitializeSmootherState();
    x_sm.assign(N, VecXd());
    P_sm.assign(N, MatXd());
    x_sm[N - 1] = filt.GetSmoothedState(N - 1);
    P_sm[N - 1] = filt.GetSmoothedCovariance(N - 1);
    for (int i = N - 2; i >= 0; --i) {
      filt.UpdateSmoother(i);
      x_sm[i] = filt.GetSmoothedState(i);
      P_sm[i] = filt.GetSmoothedCovariance(i);
    }
  }
}  // namespace

TEST_CASE("filters.srif.factorization_reconstructs_covariance") {
  // The information square root R must satisfy R^T R = P^{-1}, i.e. P = R^{-1} R^{-T}.
  MatXd P(3, 3);
  P << 4.0, 0.5, -0.2, 0.5, 2.0, 0.3, -0.2, 0.3, 1.5;

  SRIF srif;
  srif.SetState(State(VecXd::Zero(3)));
  srif.SetCovariance(P);
  MatXd R = srif.GetInfoSqrt();
  MatXd P_reconstructed = (R.transpose() * R).inverse();
  REQUIRE_THAT((P - P_reconstructed).norm(), WithinAbs(0.0, 1.0e-12));
  // R is upper triangular.
  REQUIRE_THAT(R.triangularView<Eigen::StrictlyLower>().toDenseMatrix().norm(),
               WithinAbs(0.0, 1.0e-12));
}

TEST_CASE("filters.srif.matches_ekf_for_linear_measurement") {
  // A square-root information filter is an algebraic reformulation of the Kalman filter,
  // so on a linear-Gaussian problem it must reproduce the EKF's state and covariance.
  State x0(3);
  x0 << 1.0, -2.0, 0.5;
  MatXd P0(3, 3);
  P0 << 4.0, 0.5, -0.2, 0.5, 2.0, 0.3, -0.2, 0.3, 1.5;

  MatXd F(3, 3);
  F << 1.0, 0.1, 0.0, 0.0, 1.0, 0.2, 0.0, 0.0, 1.0;  // non-trivial STM
  MatXd Q = MatXd::Zero(3, 3);

  MatXd H(2, 3);
  H << 1.0, 0.2, 0.0, -0.3, 1.0, 0.1;
  MatXd R = MatXd::Zero(2, 2);
  R(0, 0) = 0.5;
  R(1, 1) = 0.8;
  VecXd z(2);
  z << 2.0, -1.0;

  auto fd = LinearDynamics(F);
  auto fp = ConstProcessNoise(Q);
  auto fm = LinearMeasurement(H, R);

  EKF ekf;
  ekf.SetState(x0);
  ekf.SetCovariance(P0);
  ekf.SetDynamicsFunction(fd);
  ekf.SetProcessNoiseFunction(fp);
  ekf.SetMeasurementFunction(fm);
  ekf.SetOutlierThreshold(1.0e12);
  ekf.Predict(1.0);
  ekf.Update(z);

  SRIF srif;
  srif.SetState(x0);
  srif.SetCovariance(P0);
  srif.SetDynamicsFunction(fd);
  srif.SetProcessNoiseFunction(fp);
  srif.SetMeasurementFunction(fm);
  srif.SetOutlierThreshold(1.0e12);
  srif.Predict(1.0);
  srif.Update(z);

  REQUIRE_THAT((ekf.GetState() - srif.GetState()).norm().val(), WithinAbs(0.0, 1.0e-9));
  REQUIRE_THAT((ekf.GetCovariance() - srif.GetCovariance()).norm(), WithinAbs(0.0, 1.0e-9));
}

TEST_CASE("filters.srif.matches_ekf_with_process_noise") {
  // Same equivalence, now exercising the information-form time update with a dense,
  // full-rank process noise (P_predict = F P F^T + Q).
  State x0(3);
  x0 << 0.0, 0.0, 0.0;
  MatXd P0 = 3.0 * MatXd::Identity(3, 3);

  MatXd F(3, 3);
  F << 1.0, 0.2, 0.0, 0.0, 1.0, 0.1, 0.05, 0.0, 1.0;
  MatXd L(3, 3);
  L << 0.3, 0.0, 0.0, 0.1, 0.4, 0.0, 0.05, 0.2, 0.25;
  MatXd Q = L * L.transpose();  // dense SPD process noise
  REQUIRE(std::abs(Q(0, 1)) > 1.0e-6);

  MatXd H(1, 3);
  H << 1.0, 0.5, -0.2;
  MatXd R(1, 1);
  R(0, 0) = 0.4;
  VecXd z(1);
  z << 1.5;

  auto fd = LinearDynamics(F);
  auto fp = ConstProcessNoise(Q);
  auto fm = LinearMeasurement(H, R);

  EKF ekf;
  ekf.SetState(x0);
  ekf.SetCovariance(P0);
  ekf.SetDynamicsFunction(fd);
  ekf.SetProcessNoiseFunction(fp);
  ekf.SetMeasurementFunction(fm);
  ekf.SetOutlierThreshold(1.0e12);
  ekf.Predict(1.0);
  ekf.Update(z);

  SRIF srif;
  srif.SetState(x0);
  srif.SetCovariance(P0);
  srif.SetDynamicsFunction(fd);
  srif.SetProcessNoiseFunction(fp);
  srif.SetMeasurementFunction(fm);
  srif.SetOutlierThreshold(1.0e12);
  srif.Predict(1.0);
  srif.Update(z);

  REQUIRE_THAT((ekf.GetState() - srif.GetState()).norm().val(), WithinAbs(0.0, 1.0e-9));
  REQUIRE_THAT((ekf.GetCovariance() - srif.GetCovariance()).norm(), WithinAbs(0.0, 1.0e-9));
}

TEST_CASE("filters.srif.smoother_matches_ekf") {
  // The SRIF inherits the EKF's RTS smoother, and populates the prior/posterior
  // covariances and STMs it consumes identically -- so the smoothed trajectory must
  // match the EKF's across the whole arc.
  State x0(3);
  x0 << 0.5, -1.0, 0.2;
  MatXd P0 = 2.0 * MatXd::Identity(3, 3);

  MatXd F(3, 3);
  F << 1.0, 0.15, 0.0, 0.0, 1.0, 0.1, 0.0, 0.0, 1.0;
  MatXd Q = 0.01 * MatXd::Identity(3, 3);

  MatXd H(1, 3);
  H << 1.0, 0.0, 0.3;
  MatXd R(1, 1);
  R(0, 0) = 0.25;

  std::vector<double> times{0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<VecXd> zs;
  for (double v : {0.4, 0.9, 1.3, 0.7, -0.2}) {
    VecXd z(1);
    z << v;
    zs.push_back(z);
  }

  auto fd = LinearDynamics(F);
  auto fp = ConstProcessNoise(Q);
  auto fm = LinearMeasurement(H, R);

  std::vector<VecXd> x_ekf, x_srif;
  std::vector<MatXd> P_ekf, P_srif;
  EKF ekf;
  SRIF srif;
  RunAndSmooth(ekf, x0, P0, fd, fp, fm, times, zs, x_ekf, P_ekf);
  RunAndSmooth(srif, x0, P0, fd, fp, fm, times, zs, x_srif, P_srif);

  for (size_t i = 0; i < times.size(); ++i) {
    REQUIRE_THAT((x_ekf[i] - x_srif[i]).norm(), WithinAbs(0.0, 1.0e-8));
    REQUIRE_THAT((P_ekf[i] - P_srif[i]).norm(), WithinAbs(0.0, 1.0e-8));
    // Smoothed covariance must be symmetric positive-definite.
    Eigen::SelfAdjointEigenSolver<MatXd> es(P_srif[i]);
    REQUIRE(es.eigenvalues().minCoeff() > 0.0);
  }
}

TEST_CASE("filters.srif.cwna_process_noise") {
  // Continuous white-noise acceleration: a 6x6 [pos; vel] block matrix with the
  // standard van-Loan integrals q3=psd*dt^3/3, q2=psd*dt^2/2, q1=psd*dt.
  const double dt = 5.0, psd = 2.5e-3;
  MatXd Q = CwnaProcessNoise(dt, psd);

  REQUIRE(Q.rows() == 6);
  REQUIRE(Q.cols() == 6);
  const double q3 = psd * dt * dt * dt / 3.0;
  const double q2 = psd * dt * dt / 2.0;
  const double q1 = psd * dt;
  for (int k = 0; k < 3; ++k) {
    REQUIRE_THAT(Q(k, k), WithinRel(q3, 1e-12));
    REQUIRE_THAT(Q(k, k + 3), WithinRel(q2, 1e-12));
    REQUIRE_THAT(Q(k + 3, k), WithinRel(q2, 1e-12));
    REQUIRE_THAT(Q(k + 3, k + 3), WithinRel(q1, 1e-12));
  }
  // Off-block cross terms between distinct axes are zero, and Q is symmetric PSD.
  REQUIRE_THAT((Q - Q.transpose()).norm(), WithinAbs(0.0, 1e-15));
  REQUIRE_THAT(Q(0, 1), WithinAbs(0.0, 1e-15));
  Eigen::SelfAdjointEigenSolver<MatXd> es(Q);
  REQUIRE(es.eigenvalues().minCoeff() >= 0.0);
}
