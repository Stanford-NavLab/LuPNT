#include "lupnt/applications/angles_odts/angles_odts_app.h"

#include <cmath>

#include "lupnt/agents/spacecraft.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/sat_bearing_measurement.h"
#include "lupnt/numerics/filters/filter_utils.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  namespace {
    constexpr double kArcsecToRad = 4.84813681109536e-06;  // pi / (180 * 3600)
  }  // namespace

  AnglesOdtsApp::AnglesOdtsApp(Config& config) : Application(config) {
    auto g = [&](const char* k, double d) { return config[k].as<double>(d); };
    auto gi = [&](const char* k, int d) { return config[k].as<int>(d); };
    auto gb = [&](const char* k, bool d) { return config[k].as<bool>(d); };

    LUPNT_CHECK(config["target"], "AnglesOdtsApp requires a `target` agent name", "AnglesOdtsApp");
    target_name_ = config["target"].as<std::string>();

    seed_ = gi("seed", seed_);
    dt_s_ = g("dt_s", dt_s_);
    duration_s_ = g("duration_s", duration_s_);
    monte_carlo_runs_ = gi("monte_carlo_runs", monte_carlo_runs_);

    // Bearing noise: `angle_sigma_rad` takes precedence, else `angle_sigma_arcsec`.
    if (config["angle_sigma_rad"]) {
      angle_sigma_rad_ = g("angle_sigma_rad", angle_sigma_rad_);
    } else if (config["angle_sigma_arcsec"]) {
      angle_sigma_rad_ = g("angle_sigma_arcsec", 5.0) * kArcsecToRad;
    }

    // Filter force model: unified `force_model:` block (bodies list), scalar keys as fallback.
    if (config["force_model"]) {
      const ForceModelSpec fm = ParseForceModelSpec(config["force_model"]);
      moon_gravity_degree_filter_ = fm.moon_degree;
      moon_gravity_order_filter_ = fm.moon_order;
      include_earth_ = fm.include_earth;
      include_sun_ = fm.include_sun;
      use_relativity_ = fm.relativity;
    } else {
      moon_gravity_degree_filter_ = gi("moon_gravity_degree_filter", moon_gravity_degree_filter_);
      moon_gravity_order_filter_ = gi("moon_gravity_order_filter", moon_gravity_order_filter_);
      include_earth_ = gb("include_earth", include_earth_);
      include_sun_ = gb("include_sun", include_sun_);
      use_relativity_ = gb("use_relativity", use_relativity_);
    }
    integration_step_s_ = g("integration_step_s", integration_step_s_);
    process_accel_sigma_mps2_ = g("process_accel_sigma_mps2", process_accel_sigma_mps2_);
    initial_position_sigma_m_ = g("initial_position_sigma_m", initial_position_sigma_m_);
    initial_velocity_sigma_mps_ = g("initial_velocity_sigma_mps", initial_velocity_sigma_mps_);
    outlier_threshold_ = g("outlier_threshold", outlier_threshold_);
  }

  Ptr<NBodyDynamics> AnglesOdtsApp::BuildFilterDynamics() const {
    auto orbit = MakePtr<NBodyDynamics>();
    orbit->SetIntegrator(IntegratorType::RKF45);
    orbit->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
    orbit->AddBody(Body::Moon(moon_gravity_degree_filter_, moon_gravity_order_filter_));
    if (include_earth_) orbit->AddBody(Body::Earth());
    if (include_sun_) orbit->AddBody(Body::Sun());
    orbit->SetFrame(Frame::MOON_CI);
    orbit->SetTimeStep(integration_step_s_);
    orbit->SetAutodiff(true);  // the filter needs its analytic state-transition matrix
    orbit->SetUseRelativity(use_relativity_);
    return orbit;
  }

  void AnglesOdtsApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "AnglesOdtsApp");
    Simulation* sim = agent_->GetSimulation();
    epoch0_ = GetLupntEpoch();

    observer_ = dynamic_cast<Spacecraft*>(agent_);
    LUPNT_CHECK(observer_, "AnglesOdtsApp must run on a Spacecraft", "AnglesOdtsApp");
    target_ = dynamic_cast<Spacecraft*>(sim->GetAgent(target_name_));
    LUPNT_CHECK(target_, fmt::format("target `{}` is not a Spacecraft", target_name_),
                "AnglesOdtsApp");

    // Single end-of-arc solve (DEVICE priority, so it runs after the agents' truth Steps at the
    // final epoch). The whole arc is filtered here so a Monte-Carlo ensemble runs in one sim.
    sim->Schedule(
        sim->GetDuration(), [this](Real) { Solve(); }, Event::SINGLE_EVENT,
        Event::Priority::DEVICE);
  }

  void AnglesOdtsApp::Solve() {
    if (solved_) return;
    solved_ = true;

    const double duration
        = duration_s_ > 0.0 ? duration_s_ : agent_->GetSimulation()->GetDuration().val();
    const int N = static_cast<int>(std::floor(duration / dt_s_ + 1.0e-9)) + 1;
    t_grid_.resize(N);
    for (int k = 0; k < N; ++k) t_grid_(k) = k * dt_s_;

    // ---- Truth trajectories (observer + target), forward-propagated from the agents' own
    // truth dynamics captured at t=0. 8-state [r,v,cb,cd]; only [r,v] (head 6) is used. ----
    MatXd obs_truth(N, kState), tgt_pos(N, 3);
    State xo = observer_->GetTruthStateAt(0.0);
    State xt = target_->GetTruthStateAt(0.0);
    obs_truth.row(0) = xo.head(kState).cast<double>().transpose();
    tgt_pos.row(0) = xt.head(3).cast<double>().transpose();
    for (int k = 1; k < N; ++k) {
      xo = observer_->GetTruthDynamics()->Propagate(xo, t_grid_(k - 1), t_grid_(k), nullptr);
      xt = target_->GetTruthDynamics()->Propagate(xt, t_grid_(k - 1), t_grid_(k), nullptr);
      obs_truth.row(k) = xo.head(kState).cast<double>().transpose();
      tgt_pos.row(k) = xt.head(3).cast<double>().transpose();
    }
    truth_state_ = obs_truth;

    // ---- Monte-Carlo ensemble of EKF runs (independent initial-error + measurement draws). ----
    const int MC = std::max(1, monte_carlo_runs_);
    VecXd sse_pos = VecXd::Zero(N), sse_vel = VecXd::Zero(N);

    for (int run = 0; run < MC; ++run) {
      std::mt19937 rng(static_cast<unsigned int>(seed_) + 1000u * static_cast<unsigned int>(run));

      auto dyn = BuildFilterDynamics();
      EKF ekf;
      ekf.SetDynamicsFunction(
          [dyn](const State& x, Real t0, Real tf, const State*, MatXd* F) -> State {
            MatXd stm;
            State xf = dyn->Propagate(x, t0, tf, nullptr, &stm);
            if (F) *F = stm;
            return xf;
          });
      const double sa = process_accel_sigma_mps2_;
      ekf.SetProcessNoiseFunction([sa](const State& /*x*/, Real t0, Real tf) -> MatXd {
        const double dt = std::abs((tf - t0).val());
        const Mat3d Qacc = (sa * sa) * Mat3d::Identity();
        return MatXd(ProcessNoisePosVel(Qacc, dt));
      });
      ekf.SetOutlierThreshold(outlier_threshold_);

      // Seed: observer truth at t=0 + a random initial-error draw.
      VecXd x0 = obs_truth.row(0).transpose();
      for (int j = 0; j < 3; ++j) x0(j) += SampleNormal(0.0, initial_position_sigma_m_, &rng).val();
      for (int j = 3; j < 6; ++j)
        x0(j) += SampleNormal(0.0, initial_velocity_sigma_mps_, &rng).val();
      MatXd P0 = MatXd::Zero(kState, kState);
      for (int j = 0; j < 3; ++j) P0(j, j) = initial_position_sigma_m_ * initial_position_sigma_m_;
      for (int j = 3; j < 6; ++j)
        P0(j, j) = initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_;
      State x0s(kState);
      x0s = x0.cast<Real>();
      x0s.SetFrame(Frame::MOON_CI);
      ekf.SetTime(0.0);
      ekf.SetState(x0s);
      ekf.SetCovariance(P0);

      if (run == 0) {
        est_state_ = MatXd::Zero(N, kState);
        sigma_state_ = MatXd::Zero(N, kState);
        est_state_.row(0) = ekf.GetState().cast<double>().transpose();
        sigma_state_.row(0) = ekf.GetCovariance().diagonal().cwiseSqrt().transpose();
      }
      {
        Vec3d dr = ekf.GetState().cast<double>().head(3) - obs_truth.row(0).head(3).transpose();
        Vec3d dv = ekf.GetState().cast<double>().segment(3, 3)
                   - obs_truth.row(0).segment(3, 3).transpose();
        sse_pos(0) += dr.squaredNorm();
        sse_vel(0) += dv.squaredNorm();
      }

      for (int k = 1; k < N; ++k) {
        ekf.Predict(t_grid_(k));

        const Vec3d r_obs_true = obs_truth.row(k).head(3).transpose();
        const Vec3d r_tgt = tgt_pos.row(k).transpose();
        const Vec3d u_true = (r_tgt - r_obs_true).normalized();

        SatBearingMeasurement::Config mc;
        mc.idx_position = 0;
        mc.target_pos_mci = r_tgt;
        mc.sigma_rad = angle_sigma_rad_;
        SatBearingMeasurement meas(mc);
        ekf.SetMeasurementFunction(meas.CreateFunction());

        VecX z(3);
        for (int j = 0; j < 3; ++j)
          z(j) = u_true(j) + SampleNormal(0.0, angle_sigma_rad_, &rng).val();
        ekf.Update(z);

        VecXd xk = ekf.GetState().cast<double>();
        Vec3d dr = xk.head(3) - r_obs_true;
        Vec3d dv = xk.segment(3, 3) - obs_truth.row(k).segment(3, 3).transpose();
        sse_pos(k) += dr.squaredNorm();
        sse_vel(k) += dv.squaredNorm();
        if (run == 0) {
          est_state_.row(k) = xk.transpose();
          sigma_state_.row(k) = ekf.GetCovariance().diagonal().cwiseSqrt().transpose();
        }
      }
    }

    // ---- Ensemble statistics ----
    pos_err_rms_.resize(N);
    vel_err_rms_.resize(N);
    for (int k = 0; k < N; ++k) {
      pos_err_rms_(k) = std::sqrt(sse_pos(k) / MC);
      vel_err_rms_(k) = std::sqrt(sse_vel(k) / MC);
    }
    final_pos_err_m_ = pos_err_rms_(N - 1);
    // RMS over the converged tail (last quarter of the arc).
    const int k0 = std::max(1, 3 * N / 4);
    double acc = 0.0;
    int cnt = 0;
    for (int k = k0; k < N; ++k) {
      acc += pos_err_rms_(k) * pos_err_rms_(k);
      cnt++;
    }
    rms_pos_err_m_ = cnt > 0 ? std::sqrt(acc / cnt) : final_pos_err_m_;

    // ---- Convenience [N x 19] trajectory block: [t, truth(6), est(6), sigma(6)] (run 0). ----
    trajectory_ = MatXd::Zero(N, 1 + 3 * kState);
    for (int k = 0; k < N; ++k) {
      trajectory_(k, 0) = t_grid_(k);
      trajectory_.block(k, 1, 1, kState) = truth_state_.row(k);
      trajectory_.block(k, 1 + kState, 1, kState) = est_state_.row(k);
      trajectory_.block(k, 1 + 2 * kState, 1, kState) = sigma_state_.row(k);
    }

    Logger::Info(
        fmt::format("{}: angles-only ODTS over {} epochs, {} MC run(s): final position error "
                    "{:.3f} m (RMS tail {:.3f} m)",
                    name_, N, MC, final_pos_err_m_, rms_pos_err_m_),
        "AnglesOdtsApp");
  }

  REGISTER_FACTORY_CLASS(Application, AnglesOdtsApp)

}  // namespace lupnt
