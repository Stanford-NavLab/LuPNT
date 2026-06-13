#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  FilterDynamicsFunction IdentityDynamics() {
    return [](const State& x, Real, Real, const State*, MatXd* F) {
      if (F != nullptr) *F = MatXd::Identity(x.size(), x.size());
      return x;
    };
  }

  ProcessNoiseFunction ZeroProcessNoise() {
    return [](const State& x, Real, Real) { return MatXd::Zero(x.size(), x.size()); };
  }

  FilterMeasurementFunction LinearMeasurement(const MatXd& H, const MatXd& R) {
    return [H, R](const State& x, MatXd* H_out, MatXd* R_out) {
      if (H_out != nullptr) *H_out = H;
      if (R_out != nullptr) *R_out = R;
      return VecXd(H * x.cast<double>());
    };
  }
}  // namespace

TEST_CASE("filters.udu.factorization_reconstructs_covariance") {
  MatXd P(3, 3);
  P << 4.0, 0.5, -0.2, 0.5, 2.0, 0.3, -0.2, 0.3, 1.5;

  auto [D, U] = UDUDecomposition(P);
  MatXd P_reconstructed = UDUReconstruct(U, D);

  REQUIRE_THAT((P - P_reconstructed).norm(), WithinAbs(0.0, 1.0e-12));
}

TEST_CASE("filters.udu.update_matches_ekf_for_linear_measurement") {
  State x0(2);
  x0 << 1.0, -2.0;
  MatXd P0(2, 2);
  P0 << 4.0, 0.5, 0.5, 2.0;

  MatXd H(2, 2);
  H << 1.0, 0.2, -0.3, 1.0;
  MatXd R = MatXd::Zero(2, 2);
  R(0, 0) = 0.5;
  R(1, 1) = 0.8;
  VecXd z(2);
  z << 2.0, -1.0;

  EKF ekf;
  ekf.SetState(x0);
  ekf.SetCovariance(P0);
  ekf.SetDynamicsFunction(IdentityDynamics());
  ekf.SetProcessNoiseFunction(ZeroProcessNoise());
  ekf.SetMeasurementFunction(LinearMeasurement(H, R));
  ekf.Predict(1.0);
  ekf.Update(z);

  UDUEKF udu;
  udu.SetState(x0);
  udu.SetCovariance(P0);
  udu.SetDynamicsFunction(IdentityDynamics());
  udu.SetProcessNoiseFunction(ZeroProcessNoise());
  udu.SetMeasurementFunction(LinearMeasurement(H, R));
  udu.Predict(1.0);
  udu.Update(z);

  REQUIRE_THAT((ekf.GetState() - udu.GetState()).norm().val(), WithinAbs(0.0, 1.0e-10));
  REQUIRE_THAT((ekf.GetCovariance() - udu.GetCovariance()).norm(), WithinAbs(0.0, 1.0e-10));
}
