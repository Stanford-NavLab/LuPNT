#include "lupnt/simulations/isl_odts/isl_odts_simulation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "lupnt/applications/isl_odts_app.h"
#include "lupnt/applications/lunanet_sat_app.h"
#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"
#include "lupnt/numerics/filters/filter_utils.h"

namespace lupnt {

  namespace {

    // Per-satellite filter state layout: [own(8), consider(8) x (n_sat - 1)].
    constexpr int kSubStateSize = 8;

    Ptr<NBodyDynamics> MakeOrbitDynamics(int moon_degree, int moon_order, bool use_autodiff,
                                         const IslOdtsConfig& cfg) {
      auto dynamics = MakePtr<NBodyDynamics>();
      dynamics->SetIntegrator(IntegratorType::RKF45);
      dynamics->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      dynamics->AddBody(Body::Moon(moon_degree, moon_order));
      if (cfg.include_earth) dynamics->AddBody(Body::Earth());
      if (cfg.include_sun) dynamics->AddBody(Body::Sun());
      dynamics->SetFrame(Frame::MOON_CI);
      dynamics->SetTimeStep(cfg.integration_step_s);
      dynamics->SetAutodiff(use_autodiff);
      dynamics->SetUseRelativity(cfg.use_relativity);
      return dynamics;
    }

    Ptr<JointOrbitClockDynamics> MakeJointDynamics(Ptr<NumericalDynamics> orbit_dynamics,
                                                   bool add_clock_noise, int seed,
                                                   const IslOdtsConfig& cfg) {
      auto clock = MakePtr<ClockDynamics>();
      clock->SetModel(ClockModel::OCXO);
      clock->SetClockBiasUnit(ClockBiasUnit::SECONDS);
      clock->SetAddNoise(add_clock_noise);
      clock->SetSeed(seed);

      auto joint = MakePtr<JointOrbitClockDynamics>();
      joint->SetOrbitDynamics(orbit_dynamics);
      joint->SetClockDynamics(clock);
      joint->SetUseClockRelativity(true);
      joint->SetRelativityCenterBody(BodyId::MOON);
      joint->SetAddClockNoise(add_clock_noise);
      joint->SetFrame(Frame::MOON_CI);
      joint->SetTimeStep(cfg.integration_step_s);
      joint->SetIntegrator(IntegratorType::RKF45);
      joint->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      return joint;
    }

    State MakeInitialTruthState(const IslOdtsSatelliteConfig& sat) {
      Cart6 orbit0(sat.r0_m.cast<Real>(), sat.v0_mps.cast<Real>(), Frame::MOON_CI);
      ClockState2 clock0;
      clock0.b() = sat.clock_bias_s;
      clock0.d() = sat.clock_drift_sps;
      return State(JointOrbitClockState(orbit0, clock0));
    }

    VecXd SampleErrorPosVelClock(double sigma_r, double sigma_v, double sigma_b, double sigma_d) {
      VecXd e(kSubStateSize);
      e(0) = SampleNormal(0.0, sigma_r).val();
      e(1) = SampleNormal(0.0, sigma_r).val();
      e(2) = SampleNormal(0.0, sigma_r).val();
      e(3) = SampleNormal(0.0, sigma_v).val();
      e(4) = SampleNormal(0.0, sigma_v).val();
      e(5) = SampleNormal(0.0, sigma_v).val();
      e(6) = SampleNormal(0.0, sigma_b).val();
      e(7) = SampleNormal(0.0, sigma_d).val();
      return e;
    }

    double ComputeLinkBudgetCn0Dbhz(const IslLinkBudgetConfig& lb, double range_m) {
      if (!lb.enabled || !(range_m > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      const double fspl_db = 20.0 * std::log10(4.0 * PI * range_m * lb.frequency_hz / C);
      const double received_power_dbw = lb.tx_power_dbw + lb.tx_gain_dbi + lb.rx_gain_dbi - fspl_db;
      constexpr double kBoltzmannDbwHzK = -228.6;  // 10*log10(1.380649e-23 J/K)
      const double noise_density_dbw_hz
          = kBoltzmannDbwHzK + 10.0 * std::log10(lb.system_noise_temp_k);
      return received_power_dbw - noise_density_dbw_hz;
    }

  }  // namespace

  IslOdtsSimulation::IslOdtsSimulation(IslOdtsConfig config) : config_(std::move(config)) {}

  void IslOdtsSimulation::Setup() {
    LUPNT_CHECK(config_.satellites.size() >= 2,
                "IslOdtsConfig.satellites must contain at least 2 satellites", "IslOdts");
    LUPNT_CHECK(config_.duration_s > 0.0, "duration_s must be positive", "IslOdts");
    LUPNT_CHECK(config_.dt_s > 0.0, "dt_s must be positive", "IslOdts");
    LUPNT_CHECK(config_.integration_step_s > 0.0, "integration_step_s must be positive", "IslOdts");
    LUPNT_CHECK(config_.range_sigma_m > 0.0, "range_sigma_m must be positive", "IslOdts");
    LUPNT_CHECK(config_.range_rate_sigma_mps > 0.0, "range_rate_sigma_mps must be positive",
                "IslOdts");
    if (config_.surface_station.enabled) {
      LUPNT_CHECK(config_.surface_station.pseudorange_sigma_m > 0.0,
                  "surface_station.pseudorange_sigma_m must be positive when enabled", "IslOdts");
    }
    setup_complete_ = true;
  }

  void IslOdtsSimulation::Run() {
    LUPNT_CHECK(setup_complete_, "Call Setup() before Run()", "IslOdts");
    const IslOdtsConfig& cfg = config_;
    const int n_sat = static_cast<int>(cfg.satellites.size());
    const int n_links = n_sat - 1;
    const int n_state = kSubStateSize * n_sat;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    RandomEngine::SetSeed(static_cast<unsigned int>(cfg.seed));

    const Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    const int N = static_cast<int>(std::floor(cfg.duration_s / cfg.dt_s + 1.0e-9)) + 1;

    VecXd t_s(N);
    for (int k = 0; k < N; ++k) t_s(k) = k * cfg.dt_s;

    // --- Truth dynamics: independent per satellite so stochastic clock-noise draws
    // don't interleave between satellites.
    std::vector<Ptr<JointOrbitClockDynamics>> dyn_truth(n_sat);
    for (int j = 0; j < n_sat; ++j) {
      Ptr<NBodyDynamics> orbit_truth = MakeOrbitDynamics(cfg.moon_gravity_degree_truth,
                                                         cfg.moon_gravity_order_truth, false, cfg);
      dyn_truth[j] = MakeJointDynamics(orbit_truth, true, cfg.seed + 1000 * (j + 1), cfg);
    }

    std::vector<State> x_truth(n_sat);
    std::vector<VecXd> x_truth0(n_sat);
    for (int j = 0; j < n_sat; ++j) {
      x_truth[j] = MakeInitialTruthState(cfg.satellites[j]);
      x_truth0[j] = x_truth[j].cast<double>();
    }

    // --- Per-filter state ordering. Filter j stacks [own = sat j, then the other
    // satellites in ascending global index order]. `order[j][b]` = global sat index of
    // block b; `pos[j][g]` = block index of global sat g in filter j.
    std::vector<std::vector<int>> order(n_sat), pos(n_sat, std::vector<int>(n_sat, -1));
    for (int j = 0; j < n_sat; ++j) {
      order[j].push_back(j);
      for (int g = 0; g < n_sat; ++g)
        if (g != j) order[j].push_back(g);
      for (int b = 0; b < n_sat; ++b) pos[j][order[j][b]] = b;
    }

    // --- N parallel onboard Schmidt-EKF filters, one per satellite, each hosted like
    // flight software inside its own LunaNetSatApp (mirroring the single-hub original).
    std::vector<Ptr<IslOdtsApp>> apps(n_sat);
    std::vector<Ptr<LunaNetSatApp>> sat_apps(n_sat);
    for (int j = 0; j < n_sat; ++j) {
      VecXd x0j(n_state);
      MatXd P0j = MatXd::Zero(n_state, n_state);
      for (int b = 0; b < n_sat; ++b) {
        const int g = order[j][b];
        const double sr = (b == 0) ? cfg.initial_position_sigma_m : cfg.consider_position_sigma_m;
        const double sv = (b == 0) ? cfg.initial_velocity_sigma_mps : cfg.consider_velocity_sigma_mps;
        const double sb = (b == 0) ? cfg.initial_clock_bias_sigma_s : cfg.consider_clock_bias_sigma_s;
        const double sd
            = (b == 0) ? cfg.initial_clock_drift_sigma_sps : cfg.consider_clock_drift_sigma_sps;
        x0j.segment(kSubStateSize * b, kSubStateSize)
            = x_truth0[g] + SampleErrorPosVelClock(sr, sv, sb, sd);
        P0j.block(kSubStateSize * b, kSubStateSize * b, kSubStateSize, kSubStateSize)
            = InitialCovariancePosVelClock(sr, sv, sb, sd);
      }

      IslOdtsAppParams app_params;
      app_params.n_sat = n_sat;
      app_params.range_sigma_m = cfg.range_sigma_m;
      app_params.range_rate_sigma_mps = cfg.range_rate_sigma_mps;
      app_params.pseudorange_sigma_m = cfg.surface_station.pseudorange_sigma_m;
      app_params.process_accel_sigma_mps2 = cfg.process_accel_sigma_mps2;
      apps[j] = MakePtr<IslOdtsApp>(app_params);

      // Reduced-order, deterministic, autodiff-enabled filter dynamics (one per filter).
      Ptr<NBodyDynamics> orbit_filter = MakeOrbitDynamics(
          cfg.moon_gravity_degree_filter, cfg.moon_gravity_order_filter, true, cfg);
      Ptr<JointOrbitClockDynamics> dyn_filter = MakeJointDynamics(orbit_filter, false, 0, cfg);
      apps[j]->Configure(t0_tdb, x0j, P0j, dyn_filter);

      sat_apps[j] = MakePtr<LunaNetSatApp>();
      sat_apps[j]->SetName("lunanet_sat_" + cfg.satellites[j].name);
      sat_apps[j]->AddSubApp(apps[j]);
      sat_apps[j]->Setup();
    }

    // --- Lunar surface station geometry (fixed in the Moon principal-axis frame).
    Vec3 station_bf = Vec3::Zero();  // station position, Frame::MOON_PA
    if (cfg.surface_station.enabled) {
      State lla(3);
      VecX lla_v(3);
      lla_v << cfg.surface_station.latitude_deg, cfg.surface_station.longitude_deg,
          cfg.surface_station.altitude_m;
      lla = lla_v;
      State station_cart = LatLonAltToCart(lla, R_MOON, 0.0);
      station_bf = Vec3(station_cart.head(3));
    }
    int rotation_ptr = 0;  // round-robin pointer over satellites

    // --- Result buffers.
    results_ = IslOdtsResults{};
    results_.t_s = t_s;
    results_.satellite_names.resize(n_sat);
    for (int j = 0; j < n_sat; ++j) results_.satellite_names[j] = cfg.satellites[j].name;
    results_.truth_states.assign(n_sat, MatXd(N, kSubStateSize));
    results_.est.assign(n_sat, MatXd::Zero(N, n_state));
    results_.cov_diag.assign(n_sat, MatXd::Zero(N, n_state));
    results_.range_true_m.resize(N, n_links);
    results_.range_rate_true_mps.resize(N, n_links);
    results_.range_obs_m = MatXd::Constant(N, n_links, nan);
    results_.range_rate_obs_mps = MatXd::Constant(N, n_links, nan);
    results_.range_resid_m.assign(n_sat, MatXd::Constant(N, n_links, nan));
    results_.cn0_dbhz.resize(N, n_links);
    results_.served_sat_idx = VecXd::Constant(N, -1.0);
    results_.station_pr_true_m = VecXd::Constant(N, nan);
    results_.station_pr_obs_m = VecXd::Constant(N, nan);
    results_.station_pr_resid_m = VecXd::Constant(N, nan);
    results_.station_pos_mci = MatXd::Constant(N, 3, nan);

    // Symmetric pairwise crosslink truth (one physical link per unordered pair).
    auto compute_pair_truth = [&](MatXd* pr, MatXd* prr) {
      *pr = MatXd::Zero(n_sat, n_sat);
      *prr = MatXd::Zero(n_sat, n_sat);
      for (int a = 0; a < n_sat; ++a) {
        for (int b = a + 1; b < n_sat; ++b) {
          Vec2 y = RangeAndRangeRate(VecX(x_truth[a].head(3)), VecX(x_truth[b].head(3)),
                                     VecX(x_truth[a].segment(3, 3)), VecX(x_truth[b].segment(3, 3)));
          (*pr)(a, b) = (*pr)(b, a) = y(0).val();
          (*prr)(a, b) = (*prr)(b, a) = y(1).val();
        }
      }
    };

    // Record filter estimates/covariances (reordered into global sat index order so
    // est[j] block g is filter j's estimate of satellite g).
    auto record_filters = [&](int k) {
      for (int j = 0; j < n_sat; ++j) {
        State est_state = apps[j]->GetEstimate();
        VecXd est_d = est_state.cast<double>();
        VecXd cov_d = apps[j]->GetCovariance().diagonal();
        for (int g = 0; g < n_sat; ++g) {
          const int b = pos[j][g];
          results_.est[j].block(k, kSubStateSize * g, 1, kSubStateSize)
              = est_d.segment(kSubStateSize * b, kSubStateSize).transpose();
          results_.cov_diag[j].block(k, kSubStateSize * g, 1, kSubStateSize)
              = cov_d.segment(kSubStateSize * b, kSubStateSize).transpose();
        }
      }
    };

    // --- Epoch 0: record truth + initial estimates (no measurement).
    for (int j = 0; j < n_sat; ++j)
      results_.truth_states[j].row(0) = x_truth[j].cast<double>().transpose();
    {
      MatXd pr, prr;
      compute_pair_truth(&pr, &prr);
      for (int i = 0; i < n_links; ++i) {
        results_.range_true_m(0, i) = pr(0, i + 1);
        results_.range_rate_true_mps(0, i) = prr(0, i + 1);
        results_.cn0_dbhz(0, i) = ComputeLinkBudgetCn0Dbhz(cfg.link_budget, pr(0, i + 1));
      }
    }
    record_filters(0);

    double next_exchange_s = (cfg.consider_exchange_interval_s > 0.0)
                                 ? cfg.consider_exchange_interval_s
                                 : std::numeric_limits<double>::infinity();

    auto pbar = Logger::GetProgressBar(N - 1, "ISL ODTS");
    for (int k = 1; k < N; ++k) {
      const Real tkm1 = t0_tdb + t_s(k - 1);
      const Real tk = t0_tdb + t_s(k);

      // 1. Propagate all truth states.
      for (int j = 0; j < n_sat; ++j)
        x_truth[j] = dyn_truth[j]->Propagate(JointOrbitClockState(x_truth[j]), tkm1, tk, nullptr);

      // 2. Pairwise crosslink truth + one noisy observation per unordered pair.
      MatXd pr_true, prr_true;
      compute_pair_truth(&pr_true, &prr_true);
      MatXd pr_obs = MatXd::Zero(n_sat, n_sat), prr_obs = MatXd::Zero(n_sat, n_sat);
      for (int a = 0; a < n_sat; ++a) {
        for (int b = a + 1; b < n_sat; ++b) {
          const double r = pr_true(a, b) + SampleNormal(0.0, cfg.range_sigma_m).val();
          const double rr = prr_true(a, b) + SampleNormal(0.0, cfg.range_rate_sigma_mps).val();
          pr_obs(a, b) = pr_obs(b, a) = r;
          prr_obs(a, b) = prr_obs(b, a) = rr;
        }
      }

      // 3. Surface-station rotation: serve the next visible satellite (round-robin).
      int served = -1;
      Vec3d r_station_d = Vec3d::Zero();
      double station_pr_true = nan, station_pr_obs = nan;
      if (cfg.surface_station.enabled) {
        Vec6 st6;
        st6 << station_bf, Vec3::Zero();
        Vec6 st_mci = ConvertFrame(tk, st6, Frame::MOON_PA, Frame::MOON_CI);
        r_station_d = Vec3(st_mci.head(3)).cast<double>();
        results_.station_pos_mci.row(k) = r_station_d.transpose();

        for (int attempt = 0; attempt < n_sat; ++attempt) {
          const int cand = (rotation_ptr + attempt) % n_sat;
          Vec6 xt;
          xt << x_truth[cand].head(3), x_truth[cand].segment(3, 3);
          Vec6 xt_bf = ConvertFrame(tk, xt, Frame::MOON_CI, Frame::MOON_PA);
          Cart3 r_sat_bf(Vec3(xt_bf.head(3)), Frame::MOON_PA);
          Cart3 r_gs_bf(station_bf, Frame::MOON_PA);
          State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
          const double elevation_deg = (aer(1) * DEG).val();
          if (elevation_deg > cfg.surface_station.elevation_mask_deg) {
            served = cand;
            break;
          }
        }
        rotation_ptr = (served >= 0) ? (served + 1) % n_sat : (rotation_ptr + 1) % n_sat;

        if (served >= 0) {
          const Vec3d r_s = x_truth[served].cast<double>().head(3);
          const double b_s = x_truth[served].cast<double>()(6);
          station_pr_true = (r_s - r_station_d).norm() + C * b_s;
          station_pr_obs
              = station_pr_true + SampleNormal(0.0, cfg.surface_station.pseudorange_sigma_m).val();
        }
      }

      // 4. Stage each filter's crosslink rows (to its neighbors, in filter-local order),
      // plus the station pseudorange row for the served satellite's own filter, then step.
      for (int j = 0; j < n_sat; ++j) {
        IslOdtsMeasurementEpoch meas;
        meas.crosslink_range_m.resize(n_links);
        meas.crosslink_range_rate_mps.resize(n_links);
        for (int i = 0; i < n_links; ++i) {
          const int g = order[j][i + 1];
          meas.crosslink_range_m(i) = pr_obs(j, g);
          meas.crosslink_range_rate_mps(i) = prr_obs(j, g);
        }
        if (j == served) {
          meas.anchor_pos_mci = {r_station_d};
          meas.anchor_pseudorange_m.resize(1);
          meas.anchor_pseudorange_m(0) = station_pr_obs;
        }
        apps[j]->StageMeasurements(meas);
        sat_apps[j]->Step(tk);

        VecXd resid = apps[j]->GetPrefitResidual();
        for (int i = 0; i < n_links; ++i) results_.range_resid_m[j](k, i) = resid(2 * i);
        if (j == served) results_.station_pr_resid_m(k) = resid(2 * n_links);
      }

      // 5. Consider-state exchange: each filter broadcasts its posterior own estimate
      // (mean + covariance); every other filter overwrites the matching consider block
      // and drops that block's cross-covariance (fresh, independent prior).
      if (t_s(k) + 1.0e-9 >= next_exchange_s) {
        std::vector<VecXd> own_mean(n_sat);
        std::vector<MatXd> own_cov(n_sat);
        for (int i = 0; i < n_sat; ++i) {
          own_mean[i] = apps[i]->GetEstimate().cast<double>().segment(0, kSubStateSize);
          own_cov[i] = apps[i]->GetCovariance().block(0, 0, kSubStateSize, kSubStateSize);
        }
        for (int j = 0; j < n_sat; ++j) {
          State xj = apps[j]->GetEstimate();
          MatXd Pj = apps[j]->GetCovariance();
          for (int b = 1; b < n_sat; ++b) {
            const int g = order[j][b];
            xj.segment(kSubStateSize * b, kSubStateSize) = own_mean[g].cast<Real>();
            Pj.block(kSubStateSize * b, 0, kSubStateSize, n_state).setZero();
            Pj.block(0, kSubStateSize * b, n_state, kSubStateSize).setZero();
            Pj.block(kSubStateSize * b, kSubStateSize * b, kSubStateSize, kSubStateSize)
                = own_cov[g];
          }
          apps[j]->GetFilter()->SetState(xj);
          apps[j]->GetFilter()->SetCovariance(Pj);
        }
        next_exchange_s += cfg.consider_exchange_interval_s;
      }

      // 6. Record truth, estimates, crosslink-geometry + station diagnostics.
      for (int j = 0; j < n_sat; ++j)
        results_.truth_states[j].row(k) = x_truth[j].cast<double>().transpose();
      for (int i = 0; i < n_links; ++i) {
        results_.range_true_m(k, i) = pr_true(0, i + 1);
        results_.range_rate_true_mps(k, i) = prr_true(0, i + 1);
        results_.range_obs_m(k, i) = pr_obs(0, i + 1);
        results_.range_rate_obs_mps(k, i) = prr_obs(0, i + 1);
        results_.cn0_dbhz(k, i) = ComputeLinkBudgetCn0Dbhz(cfg.link_budget, pr_true(0, i + 1));
      }
      record_filters(k);
      results_.served_sat_idx(k) = served;
      results_.station_pr_true_m(k) = station_pr_true;
      results_.station_pr_obs_m(k) = station_pr_obs;

      pbar->Update(k);
    }
    pbar->Finish();
    for (int j = 0; j < n_sat; ++j) sat_apps[j]->Finish();
  }

}  // namespace lupnt
