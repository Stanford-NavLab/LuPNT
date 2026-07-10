#include "lupnt/applications/lunar_sat_odts/ground_odts_app.h"

#include <cmath>
#include <limits>
#include <tuple>

#include "lupnt/agents/isl_satellite.h"
#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  namespace {
    Ptr<JointOrbitClockDynamics> MakeGroundDynamics(int deg, int ord, bool earth, bool sun,
                                                    bool use_rel, double step) {
      auto orbit = MakePtr<NBodyDynamics>();
      orbit->SetIntegrator(IntegratorType::RKF45);
      orbit->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      orbit->AddBody(Body::Moon(deg, ord));
      if (earth) orbit->AddBody(Body::Earth());
      if (sun) orbit->AddBody(Body::Sun());
      orbit->SetFrame(Frame::MOON_CI);
      orbit->SetTimeStep(step);
      orbit->SetAutodiff(true);
      orbit->SetUseRelativity(use_rel);
      auto clock = MakePtr<ClockDynamics>();
      clock->SetModel(ClockModel::OCXO);
      clock->SetClockBiasUnit(ClockBiasUnit::SECONDS);
      clock->SetAddNoise(false);
      auto joint = MakePtr<JointOrbitClockDynamics>();
      joint->SetOrbitDynamics(orbit);
      joint->SetClockDynamics(clock);
      joint->SetUseClockRelativity(true);
      joint->SetRelativityCenterBody(BodyId::MOON);
      joint->SetAddClockNoise(false);
      joint->SetFrame(Frame::MOON_CI);
      joint->SetTimeStep(step);
      joint->SetIntegrator(IntegratorType::RKF45);
      joint->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      return joint;
    }
  }  // namespace

  GroundOdtsApp::GroundOdtsApp(Config& config) : Application(config) {
    auto g = [&](const char* k, double d) { return config[k].as<double>(d); };
    auto gi = [&](const char* k, int d) { return config[k].as<int>(d); };
    auto gb = [&](const char* k, bool d) { return config[k].as<bool>(d); };
    seed_ = gi("seed", seed_);
    dt_s_ = g("dt_s", dt_s_);
    duration_s_ = g("duration_s", duration_s_);
    moon_gravity_degree_filter_ = gi("moon_gravity_degree_filter", moon_gravity_degree_filter_);
    moon_gravity_order_filter_ = gi("moon_gravity_order_filter", moon_gravity_order_filter_);
    include_earth_ = gb("include_earth", include_earth_);
    include_sun_ = gb("include_sun", include_sun_);
    use_relativity_ = gb("use_relativity", use_relativity_);
    integration_step_s_ = g("integration_step_s", integration_step_s_);
    pseudorange_sigma_m_ = g("pseudorange_sigma_m", pseudorange_sigma_m_);
    include_station_doppler_ = gb("enable_station_doppler", include_station_doppler_);
    station_doppler_sigma_mps_ = g("station_doppler_sigma_mps", station_doppler_sigma_mps_);
    process_accel_sigma_mps2_ = g("central_process_accel_sigma_mps2", process_accel_sigma_mps2_);
    outlier_threshold_ = g("central_outlier_threshold", outlier_threshold_);
    initial_position_sigma_m_ = g("initial_position_sigma_m", initial_position_sigma_m_);
    initial_velocity_sigma_mps_ = g("initial_velocity_sigma_mps", initial_velocity_sigma_mps_);
    initial_clock_bias_sigma_s_ = g("initial_clock_bias_sigma_s", initial_clock_bias_sigma_s_);
    initial_clock_drift_sigma_sps_
        = g("initial_clock_drift_sigma_sps", initial_clock_drift_sigma_sps_);
    if (config["satellites"])
      for (const auto& s : config["satellites"]) sat_names_.push_back(s.as<std::string>());
  }

  void GroundOdtsApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "GroundOdtsApp");
    Application::Setup();
  }

  void GroundOdtsApp::Initialize() {
    Simulation* sim = agent_->GetSimulation();
    epoch0_ = GetLupntEpoch();
    rng_.seed(static_cast<unsigned int>(seed_) + 7777u);

    LUPNT_CHECK(!sat_names_.empty(), "GroundOdtsApp requires a `satellites:` list",
                "GroundOdtsApp");
    sats_.clear();
    for (const auto& sn : sat_names_) {
      auto* s = dynamic_cast<IslSatellite*>(sim->GetAgent(sn));
      LUPNT_CHECK(s, fmt::format("`{}` is not an IslSatellite", sn), "GroundOdtsApp");
      sats_.push_back(s);
    }
    n_sat_ = static_cast<int>(sats_.size());
    n_state_ = kSub * n_sat_;

    stations_bf_.clear();
    station_mask_deg_.clear();
    if (config_["stations"]) {
      for (const auto& s : config_["stations"]) {
        Config sc(s);
        State lla(3);
        VecX v(3);
        v << sc["latitude_deg"].as<Real>(), sc["longitude_deg"].as<Real>(),
            sc["altitude_m"].as<Real>(0.0);
        lla = v;
        State cart = LatLonAltToCart(lla, R_MOON, 0.0);
        stations_bf_.push_back(Vec3(cart.head(3)));
        station_mask_deg_.push_back(sc["elevation_mask_deg"].as<double>(5.0));
      }
    }

    const int N = static_cast<int>(std::floor(duration_s_ / dt_s_ + 1.0e-9)) + 1;
    t_grid_.resize(N);
    for (int k = 0; k < N; ++k) t_grid_(k) = k * dt_s_;

    std::normal_distribution<double> nd(0.0, 1.0);
    VecXd x0(n_state_);
    MatXd P0 = MatXd::Zero(n_state_, n_state_);
    for (int j = 0; j < n_sat_; ++j) {
      VecXd tr = sats_[j]->GetTruthStateAt(0.0).cast<double>();
      VecXd e(kSub);
      for (int i = 0; i < 3; ++i) e(i) = initial_position_sigma_m_ * nd(rng_);
      for (int i = 3; i < 6; ++i) e(i) = initial_velocity_sigma_mps_ * nd(rng_);
      e(6) = initial_clock_bias_sigma_s_ * nd(rng_);
      e(7) = initial_clock_drift_sigma_sps_ * nd(rng_);
      x0.segment(kSub * j, kSub) = tr + e;
      VecXd d(kSub);
      d << initial_position_sigma_m_ * initial_position_sigma_m_,
          initial_position_sigma_m_ * initial_position_sigma_m_,
          initial_position_sigma_m_ * initial_position_sigma_m_,
          initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_,
          initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_,
          initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_,
          initial_clock_bias_sigma_s_ * initial_clock_bias_sigma_s_,
          initial_clock_drift_sigma_sps_ * initial_clock_drift_sigma_sps_;
      P0.block(kSub * j, kSub * j, kSub, kSub) = MatXd(d.asDiagonal());
    }

    auto dyn
        = MakeGroundDynamics(moon_gravity_degree_filter_, moon_gravity_order_filter_,
                             include_earth_, include_sun_, use_relativity_, integration_step_s_);
    ekf_ = MakePtr<EKF>();
    State x0s(n_state_);
    x0s = x0.cast<Real>();
    x0s.SetFrame(Frame::MOON_CI);
    ekf_->SetTime(0.0);
    ekf_->SetState(x0s);
    ekf_->SetCovariance(P0);
    const int n_sat = n_sat_;
    ekf_->SetDynamicsFunction(
        [dyn, n_sat](const State& x, Real t0, Real tf, const State*, MatXd* F) -> State {
          const int n = n_sat * kSub;
          State xf(n);
          xf.SetFrame(x.GetFrame());
          if (F) F->setZero(n, n);
          for (int j = 0; j < n_sat; ++j) {
            const int off = j * kSub;
            State sub(kSub);
            sub = x.segment(off, kSub);
            sub.SetFrame(x.GetFrame());
            JointOrbitClockState xj(sub);
            if (F) {
              MatXd Fj;
              xf.segment(off, kSub) = dyn->Propagate(xj, t0, tf, nullptr, &Fj);
              F->block(off, off, kSub, kSub) = Fj;
            } else {
              xf.segment(off, kSub) = dyn->Propagate(xj, t0, tf, nullptr);
            }
          }
          return xf;
        });
    const double sigma_a = process_accel_sigma_mps2_;
    ekf_->SetProcessNoiseFunction([n_sat, sigma_a](const State& x, Real t0, Real tf) {
      const double dt = std::abs((tf - t0).val());
      MatXd Qsub = MatXd::Zero(kSub, kSub);
      Mat3d Qacc = std::pow(sigma_a, 2) * Mat3d::Identity();
      Qsub.block(0, 0, 6, 6) = ProcessNoisePosVel(Qacc, dt);
      Qsub.block(6, 6, 2, 2) = ClockDynamics::TwoStateNoise(ClockModel::OCXO, dt).cast<double>();
      MatXd Q = MatXd::Zero(x.size(), x.size());
      for (int j = 0; j < n_sat; ++j) Q.block(j * kSub, j * kSub, kSub, kSub) = Qsub;
      return Q;
    });
    ekf_->SetOutlierThreshold(outlier_threshold_);

    est_central_ = MatXd::Zero(N, n_state_);
    cov_central_full_.assign(n_sat_, MatXd::Zero(N, kSub * kSub));
    RecordEpoch(0);
  }

  void GroundOdtsApp::Step(Real t) {
    if (!initialized_) {
      Initialize();
      initialized_ = true;
    }
    int k = static_cast<int>(std::lround(t.val() / dt_s_));
    if (k < 1 || k >= static_cast<int>(t_grid_.size())) return;
    Real epoch_abs = epoch0_ + t;

    ekf_->Predict(t);

    // Gather all visible station->satellite pseudoranges (+ Doppler) this epoch.
    std::vector<std::tuple<Vec3d, Vec3d, int>> meas;  // (r_station, v_station, sat_idx)
    std::vector<double> pr_obs, dp_obs;
    for (size_t s = 0; s < stations_bf_.size(); ++s) {
      Vec6 st6;
      st6 << stations_bf_[s], Vec3::Zero();
      Vec6 st_mci = ConvertFrame(epoch_abs, st6, Frame::MOON_PA, Frame::MOON_CI);
      Vec3d rst = Vec3(st_mci.head(3)).cast<double>();
      Vec3d vst = Vec3(st_mci.tail(3)).cast<double>();
      for (int j = 0; j < n_sat_; ++j) {
        VecXd xj = sats_[j]->GetTruthStateAt(t).cast<double>();
        Vec6 xt;
        xt << xj.head(3).cast<Real>(), xj.segment(3, 3).cast<Real>();
        Vec6 xt_bf = ConvertFrame(epoch_abs, xt, Frame::MOON_CI, Frame::MOON_PA);
        Cart3 r_sat_bf(Vec3(xt_bf.head(3)), Frame::MOON_PA);
        Cart3 r_gs_bf(stations_bf_[s], Frame::MOON_PA);
        State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
        if ((aer(1) * DEG).val() <= station_mask_deg_[s]) continue;
        Vec3d rj = xj.head(3), vj = xj.segment(3, 3);
        double bj = xj(6), dj = xj(7);
        double pr = (rj - rst).norm() + C * bj + SampleNormal(0.0, pseudorange_sigma_m_).val();
        Vec3d u = (rj - rst) / (rj - rst).norm();
        double dp = u.dot(vj - vst) + C * dj + SampleNormal(0.0, station_doppler_sigma_mps_).val();
        meas.emplace_back(rst, vst, j);
        pr_obs.push_back(pr);
        dp_obs.push_back(dp);
      }
    }

    const int np = static_cast<int>(meas.size());
    if (np > 0) {
      const bool dopp = include_station_doppler_;
      const double var_pr = pseudorange_sigma_m_ * pseudorange_sigma_m_;
      const double var_dp = station_doppler_sigma_mps_ * station_doppler_sigma_mps_;
      ekf_->SetMeasurementFunction(
          [meas, dopp, var_pr, var_dp](const State& x, MatXd* H, MatXd* R) -> VecXd {
            const int nn = static_cast<int>(meas.size());
            const int mm = dopp ? 2 * nn : nn;
            const VecXd xd = x.cast<double>();
            VecXd y(mm);
            if (H) H->setZero(mm, x.size());
            if (R) {
              *R = MatXd::Zero(mm, mm);
              for (int i = 0; i < nn; ++i) {
                (*R)(i, i) = var_pr;
                if (dopp) (*R)(nn + i, nn + i) = var_dp;
              }
            }
            for (int i = 0; i < nn; ++i) {
              const Vec3d& rst = std::get<0>(meas[i]);
              const Vec3d& vst = std::get<1>(meas[i]);
              const int j = std::get<2>(meas[i]);
              const Vec3d rj = xd.segment(kSub * j, 3);
              const Vec3d vj = xd.segment(kSub * j + 3, 3);
              const double bj = xd(kSub * j + 6), dj = xd(kSub * j + 7);
              const Vec3d dr = rj - rst;
              const double rng = dr.norm();
              const Vec3d u = dr / rng;
              y(i) = rng + C * bj;
              if (H) {
                H->block(i, kSub * j, 1, 3) = u.transpose();
                (*H)(i, kSub* j + 6) = C;
              }
              if (dopp) {
                const double rdg = u.dot(vj - vst);
                y(nn + i) = rdg + C * dj;
                if (H) {
                  H->block(nn + i, kSub * j, 1, 3) = ((vj - vst) - rdg * u).transpose() / rng;
                  H->block(nn + i, kSub * j + 3, 1, 3) = u.transpose();
                  (*H)(nn + i, kSub * j + 7) = C;
                }
              }
            }
            return y;
          });
      const int m = dopp ? 2 * np : np;
      VecX cy(m);
      for (int i = 0; i < np; ++i) {
        cy(i) = pr_obs[i];
        if (dopp) cy(np + i) = dp_obs[i];
      }
      ekf_->Update(cy);
    }
    RecordEpoch(k);
  }

  void GroundOdtsApp::RecordEpoch(int k) {
    VecXd xc = ekf_->GetState().cast<double>();
    est_central_.row(k) = xc.transpose();
    MatXd Pc = ekf_->GetCovariance();
    for (int j = 0; j < n_sat_; ++j) {
      MatXd Pj = Pc.block(kSub * j, kSub * j, kSub, kSub);
      for (int r = 0; r < kSub; ++r)
        for (int c = 0; c < kSub; ++c) cov_central_full_[j](k, kSub * r + c) = Pj(r, c);
    }
  }

  REGISTER_FACTORY_CLASS(Application, GroundOdtsApp)

}  // namespace lupnt
