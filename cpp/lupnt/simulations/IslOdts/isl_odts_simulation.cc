#include "lupnt/simulations/IslOdts/isl_odts_simulation.h"

#include <cmath>
#include <limits>

#include "lupnt/filters/filter_utils.h"
#include "lupnt/filters/schmidt_ekf.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/measurements.h"

namespace lupnt {

  namespace {

    // Consider-filter state layout: [own(8), consider_1(8), ..., consider_{L}(8)]
    // for a hub linked to L = n_sat - 1 other satellites.
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

    // Extract sub-state `x.segment(offset, 8)` as a properly labeled
    // JointOrbitClockState, suitable for propagation by a JointOrbitClockDynamics.
    JointOrbitClockState ExtractSubState(const State& x, int offset) {
      State sub(kSubStateSize);
      sub = x.segment(offset, kSubStateSize);
      sub.SetFrame(x.GetFrame());
      return JointOrbitClockState(sub);
    }

    // Block-diagonal [own(8), consider_1(8), ..., consider_L(8)] dynamics: every
    // block is propagated independently (no dynamical coupling between satellites),
    // by the same (deterministic, filter-fidelity) JointOrbitClockDynamics model.
    FilterDynamicsFunction MakeCrosslinkDynamicsFunction(Ptr<JointOrbitClockDynamics> dyn,
                                                         int n_sat) {
      return [dyn, n_sat](const State& x, Real t0, Real tf, const State*, MatXd* F) -> State {
        const int n = n_sat * kSubStateSize;
        State xf(n);
        xf.SetFrame(x.GetFrame());
        if (F != nullptr) F->setZero(n, n);
        for (int j = 0; j < n_sat; ++j) {
          const int off = j * kSubStateSize;
          JointOrbitClockState xj = ExtractSubState(x, off);
          if (F != nullptr) {
            MatXd Fj;
            State xfj = dyn->Propagate(xj, t0, tf, nullptr, &Fj);
            xf.segment(off, kSubStateSize) = xfj;
            F->block(off, off, kSubStateSize, kSubStateSize) = Fj;
          } else {
            State xfj = dyn->Propagate(xj, t0, tf, nullptr);
            xf.segment(off, kSubStateSize) = xfj;
          }
        }
        return xf;
      };
    }

    MatXd ProcessNoiseSubState(double sigma_a_mps2, Real dt) {
      MatXd Q = MatXd::Zero(kSubStateSize, kSubStateSize);
      Mat3d Q_acc = std::pow(sigma_a_mps2, 2) * Mat3d::Identity();
      Q.block(0, 0, 6, 6) = ProcessNoisePosVel(Q_acc, dt);
      Q.block(6, 6, 2, 2) = ClockDynamics::TwoStateNoise(ClockModel::OCXO, dt).cast<double>();
      return Q;
    }

    ProcessNoiseFunction MakeCrosslinkProcessNoiseFunction(double sigma_a_mps2, int n_sat) {
      return [sigma_a_mps2, n_sat](const State& x, Real t0, Real tf) {
        const double dt = std::abs((tf - t0).val());
        MatXd Q_sub = ProcessNoiseSubState(sigma_a_mps2, dt);
        MatXd Q = MatXd::Zero(x.size(), x.size());
        for (int j = 0; j < n_sat; ++j) {
          Q.block(j * kSubStateSize, j * kSubStateSize, kSubStateSize, kSubStateSize) = Q_sub;
        }
        return Q;
      };
    }

    // Two-way inter-satellite range/Doppler measurement between the hub (state
    // block 0) and each of the n_links = n_sat - 1 linked satellites (state blocks
    // 1..n_links), collapsed to the geometric relative range and range-rate: a
    // coherent two-way (transponder) link is unbiased by either end's clock offset
    // to first order, so no clock terms appear here (unlike a one-way GNSS
    // pseudorange/Doppler, see `lupnt::GnssMeasurement`). Both `Range`/`RangeRate`
    // are symmetric under swapping the two satellites of a link.
    MeasurementFunction MakeCrosslinkMeasurementFunction(double sigma_range_m,
                                                          double sigma_range_rate_mps,
                                                          int n_links) {
      return [sigma_range_m, sigma_range_rate_mps, n_links](const State& x, MatXd* R) -> State {
        VecX r0 = x.head(3);
        VecX v0 = x.segment(3, 3);
        State y(2 * n_links);
        if (R != nullptr) R->setZero(2 * n_links, 2 * n_links);
        for (int i = 0; i < n_links; ++i) {
          const int off = kSubStateSize * (i + 1);
          VecX ri = x.segment(off, 3);
          VecX vi = x.segment(off + 3, 3);
          Vec2 yi = RangeAndRangeRate(r0, ri, v0, vi);
          y(2 * i) = yi(0);
          y(2 * i + 1) = yi(1);
          if (R != nullptr) {
            (*R)(2 * i, 2 * i) = sigma_range_m * sigma_range_m;
            (*R)(2 * i + 1, 2 * i + 1) = sigma_range_rate_mps * sigma_range_rate_mps;
          }
        }
        return y;
      };
    }

    Ptr<SchmidtEKF> MakeOnboardFilter(Real t0, const VecXd& x0, const MatXd& P0,
                                      Ptr<JointOrbitClockDynamics> dyn_filter,
                                      const IslOdtsConfig& cfg, int n_sat) {
      const int n_links = n_sat - 1;
      auto filter = MakePtr<SchmidtEKF>(kSubStateSize * n_links);  // trailing consider blocks
      State x0_state(kSubStateSize * n_sat);
      x0_state = x0.cast<Real>();
      x0_state.SetFrame(Frame::MOON_CI);
      filter->SetTime(t0);
      filter->SetState(x0_state);
      filter->SetCovariance(P0);
      filter->SetDynamicsFunction(MakeCrosslinkDynamicsFunction(dyn_filter, n_sat));
      filter->SetProcessNoiseFunction(
          MakeCrosslinkProcessNoiseFunction(cfg.process_accel_sigma_mps2, n_sat));
      filter->SetMeasurementFunction(GetFilterMeasurementFunction(
          MakeCrosslinkMeasurementFunction(cfg.range_sigma_m, cfg.range_rate_sigma_mps, n_links)));
      filter->SetOutlierThreshold(1.0e12);  // outlier rejection not needed for this scenario
      return filter;
    }

    double ComputeLinkBudgetCn0Dbhz(const IslLinkBudgetConfig& lb, double range_m) {
      if (!lb.enabled || !(range_m > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      const double fspl_db = 20.0 * std::log10(4.0 * PI * range_m * lb.frequency_hz / C);
      const double received_power_dbw = lb.tx_power_dbw + lb.tx_gain_dbi + lb.rx_gain_dbi - fspl_db;
      constexpr double kBoltzmannDbwHzK = -228.6;  // 10*log10(1.380649e-23 J/K)
      const double noise_density_dbw_hz = kBoltzmannDbwHzK + 10.0 * std::log10(lb.system_noise_temp_k);
      return received_power_dbw - noise_density_dbw_hz;
    }

  }  // namespace

  IslOdtsSimulation::IslOdtsSimulation(IslOdtsConfig config) : config_(std::move(config)) {}

  void IslOdtsSimulation::Setup() {
    LUPNT_CHECK(config_.satellites.size() >= 2,
                "IslOdtsConfig.satellites must contain at least 2 satellites (1 hub + >=1 linked)",
                "IslOdts");
    LUPNT_CHECK(config_.duration_s > 0.0, "duration_s must be positive", "IslOdts");
    LUPNT_CHECK(config_.dt_s > 0.0, "dt_s must be positive", "IslOdts");
    LUPNT_CHECK(config_.integration_step_s > 0.0, "integration_step_s must be positive", "IslOdts");
    LUPNT_CHECK(config_.range_sigma_m > 0.0, "range_sigma_m must be positive", "IslOdts");
    LUPNT_CHECK(config_.range_rate_sigma_mps > 0.0, "range_rate_sigma_mps must be positive",
                "IslOdts");
    setup_complete_ = true;
  }

  void IslOdtsSimulation::Run() {
    LUPNT_CHECK(setup_complete_, "Call Setup() before Run()", "IslOdts");
    const IslOdtsConfig& cfg = config_;
    const int n_sat = static_cast<int>(cfg.satellites.size());
    const int n_links = n_sat - 1;

    RandomEngine::SetSeed(static_cast<unsigned int>(cfg.seed));

    const Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    const int N = static_cast<int>(std::floor(cfg.duration_s / cfg.dt_s + 1.0e-9)) + 1;

    VecXd t_s(N);
    for (int k = 0; k < N; ++k) t_s(k) = k * cfg.dt_s;

    // Truth dynamics: independent per satellite so stochastic clock-noise draws
    // don't interleave between satellites.
    std::vector<Ptr<JointOrbitClockDynamics>> dyn_truth(n_sat);
    for (int j = 0; j < n_sat; ++j) {
      Ptr<NBodyDynamics> orbit_truth = MakeOrbitDynamics(
          cfg.moon_gravity_degree_truth, cfg.moon_gravity_order_truth, false, cfg);
      dyn_truth[j] = MakeJointDynamics(orbit_truth, true, cfg.seed + 1000 * (j + 1), cfg);
    }

    // Onboard (filter) dynamics: lower-fidelity, deterministic, autodiff-enabled;
    // shared across every state block of the single hub filter since propagation
    // is stateless given the input.
    Ptr<NBodyDynamics> orbit_filter = MakeOrbitDynamics(
        cfg.moon_gravity_degree_filter, cfg.moon_gravity_order_filter, true, cfg);
    Ptr<JointOrbitClockDynamics> dyn_filter = MakeJointDynamics(orbit_filter, false, 0, cfg);

    std::vector<State> x_truth(n_sat);
    std::vector<VecXd> x_truth0(n_sat);
    for (int j = 0; j < n_sat; ++j) {
      x_truth[j] = MakeInitialTruthState(cfg.satellites[j]);
      x_truth0[j] = x_truth[j].cast<double>();
    }

    VecXd x0(kSubStateSize * n_sat);
    x0.segment(0, kSubStateSize)
        = x_truth0[0]
          + SampleErrorPosVelClock(cfg.initial_position_sigma_m, cfg.initial_velocity_sigma_mps,
                                   cfg.initial_clock_bias_sigma_s,
                                   cfg.initial_clock_drift_sigma_sps);
    for (int i = 0; i < n_links; ++i) {
      x0.segment(kSubStateSize * (i + 1), kSubStateSize)
          = x_truth0[i + 1]
            + SampleErrorPosVelClock(cfg.consider_position_sigma_m,
                                     cfg.consider_velocity_sigma_mps,
                                     cfg.consider_clock_bias_sigma_s,
                                     cfg.consider_clock_drift_sigma_sps);
    }

    MatXd P0 = MatXd::Zero(kSubStateSize * n_sat, kSubStateSize * n_sat);
    P0.block(0, 0, kSubStateSize, kSubStateSize) = InitialCovariancePosVelClock(
        cfg.initial_position_sigma_m, cfg.initial_velocity_sigma_mps, cfg.initial_clock_bias_sigma_s,
        cfg.initial_clock_drift_sigma_sps);
    for (int i = 0; i < n_links; ++i) {
      P0.block(kSubStateSize * (i + 1), kSubStateSize * (i + 1), kSubStateSize, kSubStateSize)
          = InitialCovariancePosVelClock(cfg.consider_position_sigma_m,
                                         cfg.consider_velocity_sigma_mps,
                                         cfg.consider_clock_bias_sigma_s,
                                         cfg.consider_clock_drift_sigma_sps);
    }

    Ptr<SchmidtEKF> filter = MakeOnboardFilter(t0_tdb, x0, P0, dyn_filter, cfg, n_sat);

    results_ = IslOdtsResults{};
    results_.t_s = t_s;
    results_.satellite_names.resize(n_sat);
    for (int j = 0; j < n_sat; ++j) results_.satellite_names[j] = cfg.satellites[j].name;
    results_.truth_states.assign(n_sat, MatXd(N, kSubStateSize));
    results_.est.resize(N, kSubStateSize * n_sat);
    results_.cov_diag.resize(N, kSubStateSize * n_sat);
    results_.range_true_m.resize(N, n_links);
    results_.range_rate_true_mps.resize(N, n_links);
    results_.range_obs_m.resize(N, n_links);
    results_.range_rate_obs_mps.resize(N, n_links);
    results_.range_resid_m.resize(N, n_links);
    results_.range_rate_resid_mps.resize(N, n_links);
    results_.cn0_dbhz.resize(N, n_links);

    auto record = [&](int k, bool has_measurement, const VecXd& range_true,
                      const VecXd& range_rate_true, const VecXd& range_obs,
                      const VecXd& range_rate_obs, const VecXd& prefit_resid) {
      for (int j = 0; j < n_sat; ++j) {
        results_.truth_states[j].row(k) = x_truth[j].cast<double>().transpose();
      }
      results_.est.row(k) = filter->GetState().cast<double>().transpose();
      results_.cov_diag.row(k) = filter->GetCovariance().diagonal().transpose();

      const double nan = std::numeric_limits<double>::quiet_NaN();
      for (int i = 0; i < n_links; ++i) {
        results_.range_true_m(k, i) = range_true(i);
        results_.range_rate_true_mps(k, i) = range_rate_true(i);
        results_.range_obs_m(k, i) = has_measurement ? range_obs(i) : nan;
        results_.range_rate_obs_mps(k, i) = has_measurement ? range_rate_obs(i) : nan;
        results_.range_resid_m(k, i) = has_measurement ? prefit_resid(2 * i) : nan;
        results_.range_rate_resid_mps(k, i) = has_measurement ? prefit_resid(2 * i + 1) : nan;
        results_.cn0_dbhz(k, i) = ComputeLinkBudgetCn0Dbhz(cfg.link_budget, range_true(i));
      }
    };

    auto compute_true_links = [&](VecXd* range_true, VecXd* range_rate_true) {
      range_true->resize(n_links);
      range_rate_true->resize(n_links);
      for (int i = 0; i < n_links; ++i) {
        Vec2 y_true
            = RangeAndRangeRate(VecX(x_truth[0].head(3)), VecX(x_truth[i + 1].head(3)),
                               VecX(x_truth[0].segment(3, 3)), VecX(x_truth[i + 1].segment(3, 3)));
        (*range_true)(i) = y_true(0).val();
        (*range_rate_true)(i) = y_true(1).val();
      }
    };

    {
      VecXd range_true, range_rate_true;
      compute_true_links(&range_true, &range_rate_true);
      record(0, false, range_true, range_rate_true, VecXd(), VecXd(), VecXd());
    }

    auto pbar = Logger::GetProgressBar(N - 1, "ISL ODTS");
    for (int k = 1; k < N; ++k) {
      const Real tkm1 = t0_tdb + t_s(k - 1);
      const Real tk = t0_tdb + t_s(k);

      for (int j = 0; j < n_sat; ++j) {
        x_truth[j] = dyn_truth[j]->Propagate(JointOrbitClockState(x_truth[j]), tkm1, tk, nullptr);
      }

      filter->Predict(tk);

      VecXd range_true, range_rate_true;
      compute_true_links(&range_true, &range_rate_true);

      VecXd range_obs(n_links), range_rate_obs(n_links);
      VecXd y_obs(2 * n_links);
      for (int i = 0; i < n_links; ++i) {
        range_obs(i) = range_true(i) + SampleNormal(0.0, cfg.range_sigma_m).val();
        range_rate_obs(i) = range_rate_true(i) + SampleNormal(0.0, cfg.range_rate_sigma_mps).val();
        y_obs(2 * i) = range_obs(i);
        y_obs(2 * i + 1) = range_rate_obs(i);
      }

      filter->Update(y_obs);
      const VecXd prefit_resid = filter->GetMeasurementResidual();

      record(k, true, range_true, range_rate_true, range_obs, range_rate_obs, prefit_resid);
      pbar->Update(k);
    }
    pbar->Finish();
  }

}  // namespace lupnt
