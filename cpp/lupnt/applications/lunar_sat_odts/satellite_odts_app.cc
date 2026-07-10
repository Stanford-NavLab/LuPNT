#include "lupnt/applications/lunar_sat_odts/satellite_odts_app.h"

#include <algorithm>
#include <limits>

#include "lupnt/agents/spacecraft.h"
#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  namespace {
    constexpr int kSub = 8;

    Ptr<JointOrbitClockDynamics> MakeFilterDynamics(int deg, int ord, bool earth, bool sun,
                                                    bool use_rel, double step) {
      auto orbit = MakePtr<NBodyDynamics>();
      orbit->SetIntegrator(IntegratorType::RKF45);
      orbit->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      orbit->AddBody(Body::Moon(deg, ord));
      if (earth) orbit->AddBody(Body::Earth());
      if (sun) orbit->AddBody(Body::Sun());
      orbit->SetFrame(Frame::MOON_CI);
      orbit->SetTimeStep(step);
      orbit->SetAutodiff(true);  // filter STM
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

    VecXd SampleError8(std::mt19937& rng, double sr, double sv, double sb, double sd) {
      std::normal_distribution<double> nd(0.0, 1.0);
      VecXd e(8);
      for (int i = 0; i < 3; ++i) e(i) = sr * nd(rng);
      for (int i = 3; i < 6; ++i) e(i) = sv * nd(rng);
      e(6) = sb * nd(rng);
      e(7) = sd * nd(rng);
      return e;
    }

    MatXd InitialCov8(double sr, double sv, double sb, double sd) {
      VecXd d(8);
      d << sr * sr, sr * sr, sr * sr, sv * sv, sv * sv, sv * sv, sb * sb, sd * sd;
      return MatXd(d.asDiagonal());
    }
  }  // namespace

  SatelliteOdtsApp::SatelliteOdtsApp(Config& config) : Application(config) {
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
    range_sigma_m_ = g("range_sigma_m", range_sigma_m_);
    range_rate_sigma_mps_ = g("range_rate_sigma_mps", range_rate_sigma_mps_);
    include_time_transfer_ = gb("enable_two_way_time_transfer", include_time_transfer_);
    include_frequency_transfer_
        = gb("enable_two_way_frequency_transfer", include_frequency_transfer_);
    time_transfer_sigma_m_ = g("time_transfer_sigma_m", time_transfer_sigma_m_);
    frequency_transfer_sigma_mps_
        = g("frequency_transfer_sigma_mps", frequency_transfer_sigma_mps_);
    pseudorange_sigma_m_ = g("pseudorange_sigma_m", pseudorange_sigma_m_);
    include_station_doppler_ = gb("enable_station_doppler", include_station_doppler_);
    station_doppler_sigma_mps_ = g("station_doppler_sigma_mps", station_doppler_sigma_mps_);
    process_accel_sigma_mps2_ = g("process_accel_sigma_mps2", process_accel_sigma_mps2_);
    initial_position_sigma_m_ = g("initial_position_sigma_m", initial_position_sigma_m_);
    initial_velocity_sigma_mps_ = g("initial_velocity_sigma_mps", initial_velocity_sigma_mps_);
    initial_clock_bias_sigma_s_ = g("initial_clock_bias_sigma_s", initial_clock_bias_sigma_s_);
    initial_clock_drift_sigma_sps_
        = g("initial_clock_drift_sigma_sps", initial_clock_drift_sigma_sps_);
    consider_position_sigma_m_ = g("consider_position_sigma_m", consider_position_sigma_m_);
    consider_velocity_sigma_mps_ = g("consider_velocity_sigma_mps", consider_velocity_sigma_mps_);
    consider_clock_bias_sigma_s_ = g("consider_clock_bias_sigma_s", consider_clock_bias_sigma_s_);
    consider_clock_drift_sigma_sps_
        = g("consider_clock_drift_sigma_sps", consider_clock_drift_sigma_sps_);
    consider_exchange_interval_s_
        = g("consider_exchange_interval_s", consider_exchange_interval_s_);
    exchange_use_covariance_intersection_
        = gb("exchange_use_covariance_intersection", exchange_use_covariance_intersection_);
    exchange_ci_weight_ = g("exchange_ci_weight", exchange_ci_weight_);
    if (config["neighbors"])
      for (const auto& n : config["neighbors"]) neighbor_names_.push_back(n.as<std::string>());
  }

  void SatelliteOdtsApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "SatelliteOdtsApp");
    Simulation* sim = agent_->GetSimulation();
    self_ = dynamic_cast<Spacecraft*>(agent_);
    LUPNT_CHECK(self_, "SatelliteOdtsApp must run on an Spacecraft", "SatelliteOdtsApp");
    sat_name_ = agent_->GetName();

    // Subscribe to the constellation posterior bus (buffers neighbour broadcasts).
    sim->Subscribe("isl/posterior",
                   [this](const std::any& m) { OnPosterior(std::any_cast<IslPosteriorMsg>(m)); });

    // Schedule periodic measurement/filter Steps.
    Application::Setup();
  }

  void SatelliteOdtsApp::OnPosterior(const IslPosteriorMsg& msg) {
    if (msg.sender == sat_name_) return;  // ignore self
    for (const auto& nm : neighbor_names_)
      if (nm == msg.sender) {
        inbox_[msg.sender] = msg;
        return;
      }
  }

  void SatelliteOdtsApp::Initialize() {
    Simulation* sim = agent_->GetSimulation();
    epoch0_ = GetLupntEpoch();
    rng_.seed(static_cast<unsigned int>(seed_) + std::hash<std::string>{}(sat_name_) % 100000u);

    // Resolve neighbours (auto-discover other Spacecrafts if not listed) and stations.
    if (neighbor_names_.empty()) {
      // Not provided: nothing to auto-discover here without a registry; require explicit config.
      LUPNT_CHECK(false, "SatelliteOdtsApp requires a `neighbors:` list", "SatelliteOdtsApp");
    }
    neighbors_.clear();
    for (const auto& nm : neighbor_names_) {
      auto* nb = dynamic_cast<Spacecraft*>(sim->GetAgent(nm));
      LUPNT_CHECK(nb, fmt::format("neighbour `{}` is not an Spacecraft", nm), "SatelliteOdtsApp");
      neighbors_.push_back(nb);
    }
    // Normalize neighbour names to the resolved (simulation-prefixed) agent names, so they match
    // the `sender` field of broadcasts (which is the publisher's prefixed agent name).
    for (size_t i = 0; i < neighbors_.size(); ++i) neighbor_names_[i] = neighbors_[i]->GetName();
    n_blk_ = 1 + static_cast<int>(neighbors_.size());

    // Surface stations: fixed MOON_PA positions from their lat/lon/alt (as the coordinator did).
    stations_bf_.clear();
    station_mask_deg_.clear();
    if (config_["stations"]) {
      int idx = 0;
      for (const auto& s : config_["stations"]) {
        Config sc(s);
        State lla(3);
        VecX v(3);
        v << sc["latitude_deg"].as<Real>(), sc["longitude_deg"].as<Real>(),
            sc["altitude_m"].as<Real>(0.0);
        lla = v;
        State cart = LatLonAltToCart(lla, R_MOON, 0.0);
        std::string nm
            = sc["name"] ? sc["name"].as<std::string>() : ("station_" + std::to_string(idx));
        stations_bf_.emplace_back(Vec3(cart.head(3)), nm);
        station_mask_deg_.push_back(sc["elevation_mask_deg"].as<double>(5.0));
        idx++;
      }
    }

    const int N = static_cast<int>(std::floor(duration_s_ / dt_s_ + 1.0e-9)) + 1;
    t_grid_.resize(N);
    for (int k = 0; k < N; ++k) t_grid_(k) = k * dt_s_;

    // --- Onboard Schmidt-EKF: own block 0 + one consider block per neighbour. ---
    IslOdtsAppParams p;
    p.n_sat = n_blk_;
    p.range_sigma_m = range_sigma_m_;
    p.range_rate_sigma_mps = range_rate_sigma_mps_;
    p.pseudorange_sigma_m = pseudorange_sigma_m_;
    p.process_accel_sigma_mps2 = process_accel_sigma_mps2_;
    p.include_time_transfer = include_time_transfer_;
    p.time_transfer_sigma_m = time_transfer_sigma_m_;
    p.include_frequency_transfer = include_frequency_transfer_;
    p.frequency_transfer_sigma_mps = frequency_transfer_sigma_mps_;
    p.include_anchor_doppler = include_station_doppler_;
    p.anchor_doppler_sigma_mps = station_doppler_sigma_mps_;
    filter_ = MakePtr<IslOdtsApp>(p);

    // Initial estimate: block 0 = own truth + perturbation; blocks 1.. = neighbour truth + perturb.
    VecXd x0(kSub * n_blk_);
    MatXd P0 = MatXd::Zero(kSub * n_blk_, kSub * n_blk_);
    VecXd own_truth = self_->GetTruthStateAt(0.0).cast<double>();
    x0.segment(0, kSub)
        = own_truth
          + SampleError8(rng_, initial_position_sigma_m_, initial_velocity_sigma_mps_,
                         initial_clock_bias_sigma_s_, initial_clock_drift_sigma_sps_);
    P0.block(0, 0, kSub, kSub)
        = InitialCov8(initial_position_sigma_m_, initial_velocity_sigma_mps_,
                      initial_clock_bias_sigma_s_, initial_clock_drift_sigma_sps_);
    for (int b = 1; b < n_blk_; ++b) {
      VecXd nt = neighbors_[b - 1]->GetTruthStateAt(0.0).cast<double>();
      x0.segment(kSub * b, kSub)
          = nt
            + SampleError8(rng_, consider_position_sigma_m_, consider_velocity_sigma_mps_,
                           consider_clock_bias_sigma_s_, consider_clock_drift_sigma_sps_);
      P0.block(kSub * b, kSub * b, kSub, kSub)
          = InitialCov8(consider_position_sigma_m_, consider_velocity_sigma_mps_,
                        consider_clock_bias_sigma_s_, consider_clock_drift_sigma_sps_);
    }
    auto dyn_filter
        = MakeFilterDynamics(moon_gravity_degree_filter_, moon_gravity_order_filter_,
                             include_earth_, include_sun_, use_relativity_, integration_step_s_);
    filter_->Configure(0.0, x0, P0, dyn_filter);
    // Separate dynamics for propagating received (stale) neighbour broadcasts to the current epoch.
    exchange_dyn_
        = MakeFilterDynamics(moon_gravity_degree_filter_, moon_gravity_order_filter_,
                             include_earth_, include_sun_, use_relativity_, integration_step_s_);
    host_ = MakePtr<LunaNetSatApp>();
    host_->SetName("host_" + sat_name_);
    host_->AddSubApp(filter_);
    host_->Setup();

    truth_state_ = MatXd::Zero(N, kSub);
    own_est_ = MatXd::Zero(N, kSub);
    own_cov_diag_ = MatXd::Zero(N, kSub);
    own_cov_full_ = MatXd::Zero(N, kSub * kSub);
    next_exchange_s_ = (consider_exchange_interval_s_ > 0.0)
                           ? consider_exchange_interval_s_
                           : std::numeric_limits<double>::infinity();
    RecordEpoch(0);
  }

  void SatelliteOdtsApp::Step(Real t) {
    if (!initialized_) {
      Initialize();
      initialized_ = true;
    }
    int k = static_cast<int>(std::lround(t.val() / dt_s_));
    if (k < 1 || k >= static_cast<int>(t_grid_.size())) return;
    Simulation* sim = agent_->GetSimulation();
    const int n_links = n_blk_ - 1;
    Real epoch_abs = epoch0_ + t;

    // Own + neighbour truth at this epoch (8-state).
    VecXd xo = self_->GetTruthStateAt(t).cast<double>();

    IslOdtsMeasurementEpoch meas;
    meas.crosslink_range_m.resize(n_links);
    meas.crosslink_range_rate_mps.resize(n_links);
    if (include_time_transfer_) meas.crosslink_time_transfer_m.resize(n_links);
    if (include_frequency_transfer_) meas.crosslink_frequency_transfer_mps.resize(n_links);
    for (int i = 0; i < n_links; ++i) {
      VecXd xn = neighbors_[i]->GetTruthStateAt(t).cast<double>();
      Vec2 y = RangeAndRangeRate(VecX(xo.head(3).cast<Real>()), VecX(xn.head(3).cast<Real>()),
                                 VecX(xo.segment(3, 3).cast<Real>()),
                                 VecX(xn.segment(3, 3).cast<Real>()));
      meas.crosslink_range_m(i) = y(0).val() + SampleNormal(0.0, range_sigma_m_).val();
      meas.crosslink_range_rate_mps(i)
          = y(1).val() + SampleNormal(0.0, range_rate_sigma_mps_).val();
      if (include_time_transfer_)
        meas.crosslink_time_transfer_m(i)
            = C * (xo(6) - xn(6)) + SampleNormal(0.0, time_transfer_sigma_m_).val();
      if (include_frequency_transfer_)
        meas.crosslink_frequency_transfer_mps(i)
            = C * (xo(7) - xn(7)) + SampleNormal(0.0, frequency_transfer_sigma_mps_).val();
    }

    // Station aiding pseudoranges/Doppler for the stations that see this satellite.
    std::vector<Vec3d> a_pos, a_vel;
    std::vector<double> a_pr, a_dp;
    for (size_t s = 0; s < stations_bf_.size(); ++s) {
      const Vec3& r_bf = stations_bf_[s].first;
      Vec6 st6;
      st6 << r_bf, Vec3::Zero();
      Vec6 st_mci = ConvertFrame(epoch_abs, st6, Frame::MOON_PA, Frame::MOON_CI);
      Vec3d rst = Vec3(st_mci.head(3)).cast<double>();
      Vec3d vst = Vec3(st_mci.tail(3)).cast<double>();

      Vec6 xt;
      xt << xo.head(3).cast<Real>(), xo.segment(3, 3).cast<Real>();
      Vec6 xt_bf = ConvertFrame(epoch_abs, xt, Frame::MOON_CI, Frame::MOON_PA);
      Cart3 r_sat_bf(Vec3(xt_bf.head(3)), Frame::MOON_PA);
      Cart3 r_gs_bf(r_bf, Frame::MOON_PA);
      State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
      if ((aer(1) * DEG).val() <= station_mask_deg_[s]) continue;

      Vec3d rj = xo.head(3), vj = xo.segment(3, 3);
      double bj = xo(6), dj = xo(7);
      double pr = (rj - rst).norm() + C * bj + SampleNormal(0.0, pseudorange_sigma_m_).val();
      Vec3d u = (rj - rst) / (rj - rst).norm();
      double dp = u.dot(vj - vst) + C * dj + SampleNormal(0.0, station_doppler_sigma_mps_).val();
      a_pos.push_back(rst);
      a_vel.push_back(vst);
      a_pr.push_back(pr);
      a_dp.push_back(dp);
    }
    total_anchors_ += static_cast<int>(a_pos.size());
    if (!a_pos.empty()) {
      meas.anchor_pos_mci = a_pos;
      meas.anchor_pseudorange_m.resize(a_pos.size());
      for (size_t a = 0; a < a_pos.size(); ++a) meas.anchor_pseudorange_m(a) = a_pr[a];
      if (include_station_doppler_) {
        meas.anchor_vel_mci = a_vel;
        meas.anchor_doppler_mps.resize(a_pos.size());
        for (size_t a = 0; a < a_pos.size(); ++a) meas.anchor_doppler_mps(a) = a_dp[a];
      }
    }

    filter_->StageMeasurements(meas);
    host_->Step(t);

    // Broadcast own posterior EVERY step over the pub/sub bus (delivered as follow-on events,
    // so neighbours' inboxes carry a ~1-step-fresh broadcast -- minimal, easily-propagated
    // latency -- rather than one full exchange interval).
    {
      State xj0 = filter_->GetEstimate();
      MatXd Pj0 = filter_->GetCovariance();
      IslPosteriorMsg out;
      out.sender = sat_name_;
      out.t = t.val();
      out.mean = xj0.cast<double>().segment(0, kSub);
      out.cov = Pj0.block(0, 0, kSub, kSub);
      sim->Publish(t.val(), "isl/posterior", std::make_any<IslPosteriorMsg>(out));
    }

    // Consider-state exchange: at the exchange cadence, fuse the buffered neighbour broadcasts.
    if (t_grid_(k) + 1.0e-9 >= next_exchange_s_) {
      State xj = filter_->GetEstimate();
      MatXd Pj = filter_->GetCovariance();
      VecXd xjd = xj.cast<double>();
      const int n_state = kSub * n_blk_;

      // Fuse each neighbour block from the last broadcast we received from it, first
      // propagating that (~1-step-stale) broadcast forward to the current epoch -- exactly
      // as a receiver uses a broadcast ephemeris (its epoch is carried in the message).
      for (int b = 1; b < n_blk_; ++b) {
        auto it = inbox_.find(neighbor_names_[b - 1]);
        if (it == inbox_.end()) continue;  // no broadcast yet
        VecXd own_mean = it->second.mean;
        MatXd own_cov = it->second.cov;
        if (std::abs(it->second.t - t.val()) > 1e-6) {
          State nb0(kSub);
          nb0 = it->second.mean.cast<Real>();
          nb0.SetFrame(Frame::MOON_CI);
          MatXd Phi;
          State nbp
              = exchange_dyn_->Propagate(JointOrbitClockState(nb0), it->second.t, t, nullptr, &Phi);
          own_mean = nbp.cast<double>();
          own_cov = Phi * it->second.cov * Phi.transpose();
        }

        MatXd Pbb_best;
        VecXd mb_best;
        if (exchange_use_covariance_intersection_) {
          MatXd Pbb = Pj.block(kSub * b, kSub * b, kSub, kSub);
          VecXd mb = xjd.segment(kSub * b, kSub);
          VecXd sc = Pbb.diagonal().cwiseMax(1e-300).cwiseSqrt();
          VecXd sinv = sc.cwiseInverse();
          MatXd Sinv = sinv.asDiagonal();
          MatXd Yloc = (Sinv * Pbb * Sinv).inverse();
          MatXd Ybc = (Sinv * own_cov * Sinv).inverse();
          VecXd yloc = Yloc * sinv.cwiseProduct(mb);
          VecXd ybc = Ybc * sinv.cwiseProduct(own_mean);
          auto fuse = [&](double w, MatXd& Pn, VecXd& mn) {
            MatXd Yf = w * Yloc + (1.0 - w) * Ybc;
            VecXd yf = w * yloc + (1.0 - w) * ybc;
            MatXd Pt = Yf.inverse();
            Pt = 0.5 * (Pt + Pt.transpose());
            Pn = sc.asDiagonal() * Pt * sc.asDiagonal();
            mn = sc.cwiseProduct(Pt * yf);
          };
          if (exchange_ci_weight_ >= 0.0 && exchange_ci_weight_ <= 1.0) {
            fuse(exchange_ci_weight_, Pbb_best, mb_best);
          } else {
            double best = std::numeric_limits<double>::infinity();
            for (int gi = 1; gi <= 19; ++gi) {
              double w = static_cast<double>(gi) / 20.0;
              MatXd Pw;
              VecXd mw;
              fuse(w, Pw, mw);
              double obj = (Sinv * Pw * Sinv).trace();
              if (obj < best) {
                best = obj;
                Pbb_best = Pw;
                mb_best = mw;
              }
            }
          }
        } else {
          mb_best = own_mean;
          Pbb_best = own_cov;
        }
        xj.segment(kSub * b, kSub) = mb_best.cast<Real>();
        Pj.block(kSub * b, 0, kSub, n_state).setZero();
        Pj.block(0, kSub * b, n_state, kSub).setZero();
        Pj.block(kSub * b, kSub * b, kSub, kSub) = Pbb_best;
      }
      filter_->GetFilter()->SetState(xj);
      filter_->GetFilter()->SetCovariance(Pj);
      next_exchange_s_ += consider_exchange_interval_s_;
    }

    RecordEpoch(k);
  }

  void SatelliteOdtsApp::RecordEpoch(int k) {
    truth_state_.row(k) = self_->GetTruthStateAt(t_grid_(k)).cast<double>().transpose();
    VecXd est = filter_->GetEstimate().cast<double>();
    MatXd P = filter_->GetCovariance();
    own_est_.row(k) = est.segment(0, kSub).transpose();
    own_cov_diag_.row(k) = P.diagonal().segment(0, kSub).transpose();
    MatXd Pown = P.block(0, 0, kSub, kSub);
    for (int r = 0; r < kSub; ++r)
      for (int c = 0; c < kSub; ++c) own_cov_full_(k, kSub * r + c) = Pown(r, c);
  }

  REGISTER_FACTORY_CLASS(Application, SatelliteOdtsApp)

}  // namespace lupnt
