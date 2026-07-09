#include "lupnt/applications/ground_station_odts_app.h"

#include <limits>
#include <utility>

#include "lupnt/lupnt.h"

namespace lupnt {

  GroundStationOdtsApp::GroundStationOdtsApp(Config& config) : Application(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "GroundStationOdtsApp");

    LUPNT_CHECK(config["target"], "GroundStationOdtsApp requires a `target` agent name",
                "GroundStationOdtsApp");
    target_name_ = config["target"].as<std::string>();

    elevation_mask_deg_ = config["elevation_mask_deg"].as<double>(elevation_mask_deg_);
    use_range_ = config["use_range"].as<bool>(use_range_);
    use_range_rate_ = config["use_range_rate"].as<bool>(use_range_rate_);
    range_sigma_m_ = config["range_sigma_m"].as<double>(range_sigma_m_);
    range_rate_sigma_mps_ = config["range_rate_sigma_mps"].as<double>(range_rate_sigma_mps_);
    seed_ = config["seed"].as<int>(seed_);
    initial_position_sigma_m_
        = config["initial_position_sigma_m"].as<double>(initial_position_sigma_m_);
    initial_velocity_sigma_mps_
        = config["initial_velocity_sigma_mps"].as<double>(initial_velocity_sigma_mps_);
    batch_max_iterations_ = config["batch_max_iterations"].as<int>(batch_max_iterations_);
    batch_convergence_tol_ = config["batch_convergence_tol"].as<double>(batch_convergence_tol_);

    LUPNT_CHECK(use_range_ || use_range_rate_,
                "GroundStationOdtsApp needs at least one of use_range / use_range_rate",
                "GroundStationOdtsApp");
  }

  void GroundStationOdtsApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "GroundStationOdtsApp");
    Simulation* sim = agent_->GetSimulation();

    // Absolute TDB epoch of sim time t = 0 (dynamics/frames use t_tdb = GetLupntEpoch() + t).
    epoch0_ = GetLupntEpoch();

    // Target satellite agent + its dynamics (autodiff enabled for the STM).
    Agent* target_agent = sim->GetAgent(target_name_);
    target_ = dynamic_cast<AgentWithDynamics*>(target_agent);
    LUPNT_CHECK(target_, fmt::format("Target `{}` is not an AgentWithDynamics", target_name_),
                "GroundStationOdtsApp");
    target_dynamics_ = target_->GetDynamics();
    if (auto* nbody = dynamic_cast<NBodyDynamics*>(target_dynamics_)) nbody->SetAutodiff(true);

    // This ground station's fixed position in its body-fixed frame (velocity zero there).
    auto* gs = dynamic_cast<AgentWithDynamics*>(agent_);
    LUPNT_CHECK(gs, "GroundStationOdtsApp must run on an AgentWithDynamics (GroundStation)",
                "GroundStationOdtsApp");
    State gs_state = gs->GetState();
    station_r_ = gs_state.head(3);
    station_frame_ = gs_state.GetFrame();

    // Truth epoch state and a deliberately perturbed initial guess for the batch filter.
    x0_true_ = target_->GetStateAt(0.0).cast<double>();
    std::mt19937 guess_rng(static_cast<unsigned int>(seed_) + 1);
    Vec6d dx0;
    for (int j = 0; j < 3; ++j)
      dx0(j) = SampleNormal(0.0, initial_position_sigma_m_, &guess_rng).val();
    for (int j = 3; j < 6; ++j)
      dx0(j) = SampleNormal(0.0, initial_velocity_sigma_mps_, &guess_rng).val();
    x0_guess_ = x0_true_ + dx0;
    noise_rng_.seed(static_cast<unsigned int>(seed_));

    Logger::Info(fmt::format("{}: tracking `{}` (elevation mask {:.1f} deg)", name_, target_name_,
                             elevation_mask_deg_),
                 "GroundStationOdtsApp");

    // Schedule periodic measurement Steps (base class, using frequency_), then a single
    // end-of-arc batch Solve at a priority above APPLICATION so it runs after the last Step.
    Application::Setup();
    sim->Schedule(
        sim->GetDuration(), [this](Real) { Solve(); }, Event::SINGLE_EVENT,
        Event::Priority::DEVICE);
  }

  void GroundStationOdtsApp::Step(Real t) {
    Real epoch_abs = epoch0_ + t;

    // Target truth state (MOON_CI) and this station's state (MOON_CI) at this epoch.
    Vec6 xt = target_->GetStateAt(t);
    Vec6 st6;
    st6 << station_r_, Vec3::Zero();
    Vec6 st_mci = ConvertFrame(epoch_abs, st6, station_frame_, Frame::MOON_CI);

    // Topocentric elevation gate: convert the satellite to the station body-fixed frame.
    Vec6 xt_bf = ConvertFrame(epoch_abs, xt, Frame::MOON_CI, station_frame_);
    Cart3 r_sat_bf(Vec3(xt_bf.head(3)), station_frame_);
    Cart3 r_gs_bf(station_r_, station_frame_);
    State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
    double elevation_deg = (aer(1) * DEG).val();
    if (elevation_deg <= elevation_mask_deg_) return;

    // Frame-invariant range / range-rate w.r.t. the station's inertial state.
    Vec3d dr = (xt.head(3) - st_mci.head(3)).cast<double>();
    Vec3d dv = (xt.tail(3) - st_mci.tail(3)).cast<double>();
    double rho = dr.norm();
    double rho_dot = dr.dot(dv) / rho;

    double range = use_range_ ? rho + SampleNormal(0.0, range_sigma_m_, &noise_rng_).val()
                              : std::numeric_limits<double>::quiet_NaN();
    double range_rate = use_range_rate_
                            ? rho_dot + SampleNormal(0.0, range_rate_sigma_mps_, &noise_rng_).val()
                            : std::numeric_limits<double>::quiet_NaN();

    meas_t_.push_back(t.val());
    meas_range_.push_back(range);
    meas_range_rate_.push_back(range_rate);
  }

  void GroundStationOdtsApp::Solve() {
    if (solved_) return;
    solved_ = true;

    const int n_meas = static_cast<int>(meas_t_.size());
    LUPNT_CHECK(n_meas > 0,
                "GroundStationOdtsApp collected no visible measurements over the arc -- widen the "
                "duration or relax the elevation mask",
                "GroundStationOdtsApp");
    const int n_rows = (use_range_ ? 1 : 0) + (use_range_rate_ ? 1 : 0);

    // Station inertial (MOON_CI) states at every measurement epoch (x0-independent).
    std::vector<Vec6d> station_mci(n_meas);
    for (int k = 0; k < n_meas; ++k) {
      Real epoch_abs = epoch0_ + meas_t_[k];
      Vec6 st6;
      st6 << station_r_, Vec3::Zero();
      station_mci[k] = ConvertFrame(epoch_abs, st6, station_frame_, Frame::MOON_CI).cast<double>();
    }

    // Measurements + weights.
    const double w_range = 1.0 / (range_sigma_m_ * range_sigma_m_);
    const double w_range_rate = 1.0 / (range_rate_sigma_mps_ * range_rate_sigma_mps_);
    std::vector<VecXd> measurements(n_meas), weights(n_meas);
    for (int k = 0; k < n_meas; ++k) {
      VecXd y(n_rows), w(n_rows);
      int row = 0;
      if (use_range_) {
        y(row) = meas_range_[k];
        w(row) = w_range;
        row++;
      }
      if (use_range_rate_) {
        y(row) = meas_range_rate_[k];
        w(row) = w_range_rate;
        row++;
      }
      measurements[k] = y;
      weights[k] = w;
    }

    VecXd init_cov_diag(6);
    init_cov_diag.head(3).setConstant(initial_position_sigma_m_ * initial_position_sigma_m_);
    init_cov_diag.tail(3).setConstant(initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_);

    // Analytic design matrix: propagate x0 across the measurement epochs accumulating the
    // autodiff STM Phi(t_k, t0), and chain it with the closed-form range/range-rate partials
    // (evaluated in MOON_CI against the station's inertial state). No finite differencing.
    struct Cache {
      MatX6 grid_mci;              // [n_meas x 6]
      std::vector<Mat6d> stm_cum;  // [n_meas] of Phi(t_k, t0)
    };
    Cache cache;
    Dynamics* dyn = target_dynamics_;
    const std::vector<double>& meas_t = meas_t_;

    auto propagate = [&](const VecXd& x0) {
      cache.grid_mci.resize(n_meas, 6);
      cache.stm_cum.assign(n_meas, Mat6d::Identity());
      State x = Cart6(x0.cast<Real>(), Frame::MOON_CI);
      Real t_prev = 0.0;
      Mat6d phi = Mat6d::Identity();
      for (int k = 0; k < n_meas; ++k) {
        Real t_k = meas_t[k];
        if (abs(t_k - t_prev) > EPS) {
          MatXd stm_seg;
          x = dyn->Propagate(x, t_prev, t_k, nullptr, &stm_seg);
          phi = stm_seg * phi;
        }
        cache.grid_mci.row(k) = x.transpose();
        cache.stm_cum[k] = phi;
        t_prev = t_k;
      }
    };

    // The instantaneous range/range-rate geometry + design matrix live in
    // `GroundStationRangeMeasurement`; here we only supply the target/station states and
    // chain the closed-form partials with the accumulated STM `Phi(t_k, t0)`.
    MeasurementModelFunction model = [&](const VecXd& x0, int meas_idx) -> std::pair<VecXd, MatXd> {
      if (meas_idx == 0) propagate(x0);
      State xs = Cart6(cache.grid_mci.row(meas_idx).transpose(), Frame::MOON_CI);

      GroundStationRangeMeasurement::Config mcfg;
      mcfg.use_range = use_range_;
      mcfg.use_range_rate = use_range_rate_;
      mcfg.reference_state = station_mci[meas_idx];
      GroundStationRangeMeasurement meas(mcfg);

      MatXd h_obs;
      VecXd y = meas.Compute(xs, &h_obs).value;
      return {y, h_obs * cache.stm_cum[meas_idx]};
    };

    BatchFilterConfig bf_config;
    bf_config.use_weights = true;
    bf_config.use_initialization = false;
    bf_config.convergence_tol = batch_convergence_tol_;
    bf_config.max_iterations = batch_max_iterations_;

    BatchFilterResults res = RunBatchFilter(x0_guess_, init_cov_diag, measurements, weights, model,
                                            bf_config, x0_true_);

    x0_est_ = res.state_estimate;
    covariance_ = res.state_covariance;
    converged_ = res.converged;
    num_iterations_ = res.iterations;

    VecXd sigma = covariance_.diagonal().cwiseSqrt();
    double pos_err = (x0_est_.head(3) - x0_true_.head(3)).norm();
    double vel_err = (x0_est_.tail(3) - x0_true_.tail(3)).norm();
    Logger::Info(fmt::format("{}: batch OD over {} measurements, converged={} ({} iters)\n"
                             "  final position error {:.3f} m (formal 1-sigma {:.3f} m)\n"
                             "  final velocity error {:.6f} m/s (formal 1-sigma {:.6f} m/s)",
                             name_, n_meas, converged_, num_iterations_, pos_err,
                             sigma.head(3).norm(), vel_err, sigma.tail(3).norm()),
                 "GroundStationOdtsApp");

    DataLogger::Log(fmt::format("{}/x0_true", name_), x0_true_);
    DataLogger::Log(fmt::format("{}/x0_estimated", name_), x0_est_);
    DataLogger::Log(fmt::format("{}/x0_sigma", name_), sigma);
    Log(agent_->GetSimulation()->GetDuration());
  }

  void GroundStationOdtsApp::Log(Real t) {
    (void)t;
    if (!solved_) return;
    DataLogger::Log(fmt::format("{}/converged", name_), converged_ ? 1.0 : 0.0);
    DataLogger::Log(fmt::format("{}/num_measurements", name_), static_cast<double>(meas_t_.size()));
  }

  REGISTER_FACTORY_CLASS(Application, GroundStationOdtsApp)

}  // namespace lupnt
