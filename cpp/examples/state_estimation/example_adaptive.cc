// Compares fixed and adaptive process-noise strategies on a one-dimensional
// range/range-rate filtering problem.
#include <lupnt/lupnt.h>

#include "lupnt/filters/filter_utils.h"

using namespace lupnt;

int main() {
  // Output
  std::filesystem::path output_path = GetOutputDir("state_estimation") / "example_adaptive.h5";
  H5Easy::File output_file = GetH5File(output_path, true);

  // Time
  double dt = 0.1;  // [s] Step
  double tf = 240;  // [s] Final
  const int N_t = tf / dt + 1.5;
  VecX ts = Arange(0.0, tf + EPS, dt);  // [s] Span

  // Noise
  double sigma_rho = 2;          // [m] Range
  double sigma_rhodot = 0.1;     // [m/s] Range-rate
  double sigma_ap1 = sqrt(0.5);  // [m/s^2] Acceleration
  double sigma_a = 2.0;          // [m/s^2] Acceleration

  // Initial uncertainty
  double sigma0_r = 1.8;   // [m] Position
  double sigma0_v = 0.15;  // [m/s] Velocity
  double sigma0_a = 0.01;  // [m/s^2] Acceleration

  // Accelerations
  auto ap1 = [sigma_ap1](double t) { return SampleNormal(0.0, sigma_ap1); };
  auto ap2 = [](double t) { return cos(PI / 5.0 * t); };
  auto ap = ap2;

  MatXd Q_a = sigma_a * sigma_a * MatXd::Identity(1, 1);

  int N_wait = 0 * N_t;
  int N_steps = 30;

  Vec1d beta = Vec1d::Ones(1) * 1e-2;

  ProcessNoise process_noise;
  process_noise.SetNoiseLimits(1e-8, 1e8);
  process_noise.SetTimeStep(dt);
  process_noise.SetWindowSize(N_steps);
  process_noise.SetProcessNoise(Q_a);
  process_noise.SetWaitSteps(N_wait);
  process_noise.SetAlpha(2e-2);
  process_noise.SetBeta(beta);

  // Measurement function
  MeasurementFunction f_meas
      = [sigma_rho, sigma_rhodot](const State& x, MatXd* R = nullptr) -> State {
    Vec2 rv1{x(0), x(1)}, rv2{0.0, 0.0};
    Vec2 y = RangeAndRangeRate(rv1, rv2);

    if (R != nullptr) {
      *R = Vec2{sigma_rho * sigma_rho, sigma_rhodot * sigma_rhodot}.asDiagonal();
    }
    LUPNT_CHECK(!y.hasNaN(), "Measurement has NaN", "EKF");
    return y;
  };

  std::vector<ProcessNoiseAlgorithm> alg_list
      = {ProcessNoiseAlgorithm::SNC, ProcessNoiseAlgorithm::ASNC, ProcessNoiseAlgorithm::ADMC};

  for (const auto& alg : alg_list) {
    std::string alg_str = std::string(enum_name(alg));
    Logger::Info("Processing " + alg_str, "Main");

    process_noise.Reset();
    process_noise.SetAlgorithm(alg);

    int N_ekf = alg == ProcessNoiseAlgorithm::ADMC ? 3 : 2;
    int N_meas = 2;

    MatXd x_true = MatXd::Zero(N_t, 2);
    MatXd x_est = MatXd::Zero(N_t, N_ekf);
    MatXd sigma_est = MatXd::Zero(N_t, N_ekf);

    MatXd Phi;
    if (N_ekf == 2) {
      Phi = StateTransitionMatrixPosVel(dt, 1);
    } else {
      Phi = StateTransitionMatrixPosVelAcc(dt, beta);
    }

    MatXd P_est = MatXd::Zero(N_ekf, N_ekf);
    P_est(0, 0) = sigma0_r * sigma0_r;
    P_est(1, 1) = sigma0_v * sigma0_v;
    if (N_ekf == 3) {
      P_est(2, 2) = sigma0_a * sigma0_a;
    }

    RandomEngine::SetSeed(1234);
    x_true.row(0) = Vec2{0.0, 0.0};
    x_est.row(0).head(2) = SampleMvNormal(x_true.row(0), P_est, 1).transpose();
    sigma_est.row(0) = P_est.diagonal().array().sqrt();

    // Dynamics function
    DynamicsFunction f_dyn
        = [dt, Phi](const State& x, Real t0, Real tf, const State* u = nullptr) -> State {
      (void)u;
      State xf = Phi * x;
      LUPNT_CHECK(!xf.hasNaN(), "State has NaN", "EKF");
      return xf;
    };

    // Process noise function
    ProcessNoiseFunction f_proc = [&process_noise, &dt, &Q_a](const State& x, Real t0, Real tf) {
      return process_noise.ComputeProcessNoise();
    };

    // Filter
    Ptr<EKF> filter = MakePtr<EKF>();
    filter->SetDynamicsFunction(GetFilterDynamicsFunction(f_dyn));
    filter->SetMeasurementFunction(GetFilterMeasurementFunction(f_meas));
    filter->SetProcessNoiseFunction(f_proc);
    filter->SetState(x_est.row(0));
    filter->SetCovariance(P_est);

    // True state
    for (int t = 1; t < N_t; t++) {
      Real r = x_true.row(t - 1)(0);
      Real v = x_true.row(t - 1)(1);
      x_true.row(t)(0) = r + dt * v + 0.5 * dt * dt * ap(ts(t - 1));
      x_true.row(t)(1) = v + dt * ap(ts(t - 1));
    }

    auto pbar = Logger::GetProgressBar(N_t, "Filtering");
    for (int t = 1; t < N_t; t++) {
      // Predict
      filter->Predict(ts(t));

      // Update
      MatXd R_true;
      VecX y_true = f_meas(x_true.row(t), &R_true);
      y_true = SampleMvNormal(y_true, R_true).transpose();
      filter->Update(y_true);

      // Store
      x_est.row(t) = filter->GetState().transpose();
      sigma_est.row(t) = filter->GetCovariance().diagonal().array().sqrt();

      // Adaptive process noise
      process_noise.Update(filter->GetStateCorrection(), filter->GetStateCorrectionCov(),
                           filter->GetCovarianceBar(), filter->GetCovariancePost());

      pbar->Update(t);
    }
    pbar->Finish();

    Dump(output_file, fmt::format("{}/ts", alg_str), ts);
    Dump(output_file, fmt::format("{}/x_true", alg_str), x_true);
    Dump(output_file, fmt::format("{}/x_est", alg_str), x_est);
    Dump(output_file, fmt::format("{}/sigma_est", alg_str), sigma_est);
    Dump(output_file, fmt::format("{}/N_wait", alg_str), N_wait);
    Dump(output_file, fmt::format("{}/N_steps", alg_str), N_steps);
  }
}
