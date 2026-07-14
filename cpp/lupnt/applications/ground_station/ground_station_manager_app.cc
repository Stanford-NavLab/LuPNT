#include "lupnt/applications/ground_station/ground_station_manager_app.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <utility>

#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/lupnt.h"
#include "lupnt/numerics/filters/srif.h"
#include "lupnt/simulations/world.h"

namespace lupnt {

  namespace {
    void Flatten6x6(const MatXd& P, MatXd& dst, int i) {
      for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c) dst(i, r * 6 + c) = P(r, c);
    }
  }  // namespace

  GroundStationManagerApp::GroundStationManagerApp(Config& config) : Application(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "GroundStationManagerApp");

    LUPNT_CHECK(config["target"], "GroundStationManagerApp requires a `target` agent name",
                "GroundStationManagerApp");
    target_name_ = config["target"].as<std::string>();

    seed_ = config["seed"].as<int>(seed_);
    initial_position_sigma_m_
        = config["initial_position_sigma_m"].as<double>(initial_position_sigma_m_);
    initial_velocity_sigma_mps_
        = config["initial_velocity_sigma_mps"].as<double>(initial_velocity_sigma_mps_);
    batch_max_iterations_ = config["batch_max_iterations"].as<int>(batch_max_iterations_);
    batch_convergence_tol_ = config["batch_convergence_tol"].as<double>(batch_convergence_tol_);
    obs_interval_s_ = config["obs_interval_s"].as<double>(obs_interval_s_);
    run_srif_ = config["run_srif"].as<bool>(run_srif_);
    srif_use_process_noise_ = config["srif_use_process_noise"].as<bool>(srif_use_process_noise_);
    srif_accel_psd_ = config["srif_accel_psd"].as<double>(srif_accel_psd_);

    // Estimation-side correction modelling policy.
    model_shapiro_ = config["model_shapiro"].as<bool>(model_shapiro_);
    model_solid_earth_tide_ = config["model_solid_earth_tide"].as<bool>(model_solid_earth_tide_);
    tropo_cancel_fraction_ = std::clamp(
        config["troposphere_cancel_fraction"].as<double>(tropo_cancel_fraction_), 0.0, 1.0);
    iono_cancel_fraction_ = std::clamp(
        config["ionosphere_cancel_fraction"].as<double>(iono_cancel_fraction_), 0.0, 1.0);
    residual_delay_noise_scale_
        = config["residual_delay_noise_scale"].as<double>(residual_delay_noise_scale_);

    // Optional per-filter dynamics overrides (an `NBodyDynamics`-style force-model block). When
    // absent, the filters inherit the shared `world.force_model` (truth == filter).
    if (config["filter_dynamics"]) {
      filter_dyn_node_ = config["filter_dynamics"];
      has_filter_dyn_ = true;
    }
    if (config["batch_dynamics"]) {
      batch_dyn_node_ = config["batch_dynamics"];
      has_batch_dyn_ = true;
    }
    if (config["sequential_dynamics"]) {
      srif_dyn_node_ = config["sequential_dynamics"];
      has_srif_dyn_ = true;
    }
  }

  int GroundStationManagerApp::RegisterStation(const std::string& station_name) {
    station_names_.push_back(station_name);
    return static_cast<int>(station_names_.size()) - 1;
  }

  void GroundStationManagerApp::AddMeasurement(const StationMeasurement& m) { meas_.push_back(m); }

  Vec6d GroundStationManagerApp::ModeledStation(const StationMeasurement& m) const {
    Vec6d st = m.station_mci;
    if (model_solid_earth_tide_) st.head(3) += m.tide_disp_m;
    return st;
  }

  double GroundStationManagerApp::ModeledRangeCorrection(const StationMeasurement& m) const {
    double corr = tropo_cancel_fraction_ * m.tropo_delay_m + iono_cancel_fraction_ * m.iono_delay_m;
    if (model_shapiro_) corr += m.shapiro_delay_m;
    return corr;
  }

  double GroundStationManagerApp::RangeSigmaEffective(const StationMeasurement& m) const {
    double uncancelled = (1.0 - tropo_cancel_fraction_) * m.tropo_delay_m
                         + (1.0 - iono_cancel_fraction_) * m.iono_delay_m;
    double extra = residual_delay_noise_scale_ * uncancelled;
    return std::sqrt(m.range_sigma * m.range_sigma + extra * extra);
  }

  int GroundStationManagerApp::EpochIndex(double t) const {
    if (obs_interval_s_ <= 0.0) return 0;
    int i = static_cast<int>(std::llround(t / obs_interval_s_));
    if (i < 0) i = 0;
    if (i >= static_cast<int>(t_grid_.size())) i = static_cast<int>(t_grid_.size()) - 1;
    return i;
  }

  void GroundStationManagerApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "GroundStationManagerApp");
    Simulation* sim = agent_->GetSimulation();
    World* world = agent_->GetWorld();
    LUPNT_CHECK(world, "GroundStationManagerApp requires a World (define a `world:` block)",
                "GroundStationManagerApp");

    epoch0_ = GetLupntEpoch();
    dynamics_ = world->MakeDynamics();  // shared world force model (also supplies GetFrame)

    // Resolve each estimator's dynamics: an explicit per-filter override, else the shared
    // `filter_dynamics`, else the world force model (unchanged default -> truth == filter).
    auto build = [&](const Config& node) -> Ptr<NBodyDynamics> {
      Config cfg = node;
      auto dyn = MakePtr<NBodyDynamics>(cfg);
      dyn->SetFrame(world->GetFrame());
      dyn->SetAutodiff(true);  // the filter needs its analytic state-transition matrix
      return dyn;
    };
    batch_dynamics_ = has_batch_dyn_    ? build(batch_dyn_node_)
                      : has_filter_dyn_ ? build(filter_dyn_node_)
                                        : dynamics_;
    srif_dynamics_ = has_srif_dyn_     ? build(srif_dyn_node_)
                     : has_filter_dyn_ ? build(filter_dyn_node_)
                                       : dynamics_;

    // Uniform epoch grid over the whole arc, [0, duration] at obs_interval_s.
    const double duration = sim->GetDuration().val();
    const int n_epochs = static_cast<int>(std::floor(duration / obs_interval_s_)) + 1;
    t_grid_.resize(n_epochs);
    for (int i = 0; i < n_epochs; ++i) t_grid_(i) = static_cast<double>(i) * obs_interval_s_;

    // Truth epoch state, captured now (before any agent has stepped, so the target's
    // GetStateAt(0) is its initial state and needs no back-propagation).
    x0_true_ = world->GetStateAt(target_name_, 0.0).cast<double>();

    // Deliberately perturbed a-priori guess for the batch filter.
    std::mt19937 guess_rng(static_cast<unsigned int>(seed_) + 1);
    Vec6d dx0;
    for (int j = 0; j < 3; ++j)
      dx0(j) = SampleNormal(0.0, initial_position_sigma_m_, &guess_rng).val();
    for (int j = 3; j < 6; ++j)
      dx0(j) = SampleNormal(0.0, initial_velocity_sigma_mps_, &guess_rng).val();
    x0_guess_ = x0_true_ + dx0;

    // Single end-of-arc solve, at DEVICE priority so it runs after every station's
    // APPLICATION-priority tracking Step at the final epoch.
    sim->Schedule(
        sim->GetDuration(), [this](Real) { Solve(); }, Event::SINGLE_EVENT,
        Event::Priority::DEVICE);
  }

  void GroundStationManagerApp::Solve() {
    if (solved_) return;
    solved_ = true;

    const int n_meas = static_cast<int>(meas_.size());
    LUPNT_CHECK(n_meas > 0,
                "GroundStationManagerApp aggregated no measurements over the arc -- widen the "
                "duration or relax the elevation masks",
                "GroundStationManagerApp");
    const int n_epochs = static_cast<int>(t_grid_.size());

    Logger::Info(fmt::format("{}: centralized OD over {} measurements from {} stations", name_,
                             n_meas, station_names_.size()),
                 "GroundStationManagerApp");

    // Propagate an epoch state x0 over the uniform grid, returning the world-frame
    // trajectory and the cumulative STMs Phi(t_i, t0) (autodiff, chained per segment).
    auto PropagateGridStm = [&](const VecXd& x0, MatX6& grid, std::vector<Mat6d>& stm_cum) {
      grid.resize(n_epochs, 6);
      stm_cum.assign(n_epochs, Mat6d::Identity());
      State x = Cart6(x0.cast<Real>(), batch_dynamics_->GetFrame());
      grid.row(0) = x.transpose();
      Mat6d phi = Mat6d::Identity();
      for (int i = 1; i < n_epochs; ++i) {
        MatXd stm_seg;
        x = batch_dynamics_->Propagate(x, t_grid_(i - 1), t_grid_(i), nullptr, &stm_seg);
        phi = stm_seg * phi;
        stm_cum[i] = phi;
        grid.row(i) = x.transpose();
      }
    };
    auto PropagateGrid = [&](const VecXd& x0) -> MatX6 {
      MatX6 grid(n_epochs, 6);
      State x = Cart6(x0.cast<Real>(), batch_dynamics_->GetFrame());
      grid.row(0) = x.transpose();
      for (int i = 1; i < n_epochs; ++i) {
        x = batch_dynamics_->Propagate(x, t_grid_(i - 1), t_grid_(i), nullptr);
        grid.row(i) = x.transpose();
      }
      return grid;
    };

    // Closed-form range/range-rate observation for one measurement against a satellite
    // world-frame state; fills the observation partials (w.r.t. the satellite state) if h_obs.
    auto Observe = [this](const StationMeasurement& m, const Vec6d& xs, MatXd* h_obs) -> VecXd {
      const int nr = (m.has_range ? 1 : 0) + (m.has_range_rate ? 1 : 0);
      // Model the deterministic corrections into the prediction: the tide-displaced station
      // for the geometry, and the Shapiro + calibrated media delays added to the range.
      const Vec6d st = ModeledStation(m);
      Vec3d dr = xs.head(3) - st.head(3);
      Vec3d dv = xs.tail(3) - st.tail(3);
      double rho = dr.norm();
      Vec3d u = dr / rho;
      double rho_dot = dr.dot(dv) / rho;
      VecXd y(nr);
      if (h_obs) *h_obs = MatXd::Zero(nr, 6);
      int row = 0;
      if (m.has_range) {
        y(row) = rho + ModeledRangeCorrection(m);
        if (h_obs) h_obs->block(row, 0, 1, 3) = u.transpose();
        row++;
      }
      if (m.has_range_rate) {
        y(row) = rho_dot;
        if (h_obs) {
          h_obs->block(row, 0, 1, 3) = ((dv - rho_dot * u) / rho).transpose();
          h_obs->block(row, 3, 1, 3) = u.transpose();
        }
        row++;
      }
      return y;
    };

    // ---- Batch weighted least squares over all stations ----
    std::vector<VecXd> measurements(n_meas), weights(n_meas);
    for (int k = 0; k < n_meas; ++k) {
      const StationMeasurement& m = meas_[k];
      const int nr = (m.has_range ? 1 : 0) + (m.has_range_rate ? 1 : 0);
      VecXd y(nr), w(nr);
      int row = 0;
      if (m.has_range) {
        double s = RangeSigmaEffective(m);
        y(row) = m.range;
        w(row) = 1.0 / (s * s);
        row++;
      }
      if (m.has_range_rate) {
        y(row) = m.range_rate;
        w(row) = 1.0 / (m.range_rate_sigma * m.range_rate_sigma);
        row++;
      }
      measurements[k] = y;
      weights[k] = w;
    }

    MatX6 grid_mci;
    std::vector<Mat6d> stm_cum;
    MeasurementModelFunction model = [&](const VecXd& x0, int meas_idx) -> std::pair<VecXd, MatXd> {
      if (meas_idx == 0) PropagateGridStm(x0, grid_mci, stm_cum);
      const StationMeasurement& m = meas_[meas_idx];
      const int i = m.epoch_index;
      Vec6d xs = grid_mci.row(i).cast<double>().transpose();
      MatXd h_obs;
      VecXd y = Observe(m, xs, &h_obs);
      return {y, h_obs * stm_cum[i]};
    };

    VecXd init_cov_diag(6);
    init_cov_diag.head(3).setConstant(initial_position_sigma_m_ * initial_position_sigma_m_);
    init_cov_diag.tail(3).setConstant(initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_);

    BatchFilterConfig bf_config;
    bf_config.use_weights = true;
    bf_config.use_initialization = false;
    bf_config.convergence_tol = batch_convergence_tol_;
    bf_config.max_iterations = batch_max_iterations_;

    BatchFilterResults bf = RunBatchFilter(x0_guess_, init_cov_diag, measurements, weights, model,
                                           bf_config, x0_true_);
    x0_est_ = bf.state_estimate;
    covariance_ = bf.state_covariance;
    converged_ = bf.converged;
    num_iterations_ = bf.iterations;

    // ---- Iteration diagnostics ----
    const int n_iter = static_cast<int>(bf.iteration_history.size());
    iter_state_.resize(n_iter, 6);
    iter_pos_err_.resize(n_iter);
    iter_vel_err_.resize(n_iter);
    iter_corr_norm_.resize(n_iter);
    iter_weighted_rms_.resize(n_iter);
    iter_rms_range_.resize(n_iter);
    iter_rms_range_rate_.resize(n_iter);
    bool any_range = false, any_rate = false;
    for (const auto& m : meas_) {
      any_range = any_range || m.has_range;
      any_rate = any_rate || m.has_range_rate;
    }
    for (int k = 0; k < n_iter; ++k) {
      const auto& info = bf.iteration_history[k];
      iter_state_.row(k) = info.state_estimate.transpose();
      iter_corr_norm_(k) = info.correction_norm;
      iter_weighted_rms_(k) = info.weighted_rms;
      iter_pos_err_(k) = (info.state_estimate.head(3) - x0_true_.head(3)).norm();
      iter_vel_err_(k) = (info.state_estimate.tail(3) - x0_true_.tail(3)).norm();

      MatX6 grid = PropagateGrid(info.state_estimate);
      double sse_range = 0.0, sse_rate = 0.0;
      for (const auto& m : meas_) {
        Vec6d xs = grid.row(m.epoch_index).cast<double>().transpose();
        VecXd yp = Observe(m, xs, nullptr);
        int row = 0;
        if (m.has_range) sse_range += std::pow(m.range - yp(row++), 2);
        if (m.has_range_rate) sse_rate += std::pow(m.range_rate - yp(row++), 2);
      }
      iter_rms_range_(k)
          = any_range ? std::sqrt(sse_range / n_meas) : std::numeric_limits<double>::quiet_NaN();
      iter_rms_range_rate_(k)
          = any_rate ? std::sqrt(sse_rate / n_meas) : std::numeric_limits<double>::quiet_NaN();
    }

    // ---- Truth + estimated trajectories over the grid, with formal covariance ----
    truth_state_ = PropagateGrid(x0_true_).cast<double>();
    {
      MatX6 est_mci;
      std::vector<Mat6d> est_stm;
      PropagateGridStm(x0_est_, est_mci, est_stm);
      estimated_state_ = est_mci.cast<double>();
      estimated_covariance_.resize(n_epochs, 36);
      for (int i = 0; i < n_epochs; ++i) {
        Mat6d Pi = est_stm[i] * covariance_ * est_stm[i].transpose();
        Flatten6x6(Pi, estimated_covariance_, i);
      }
    }

    VecXd sigma = covariance_.diagonal().cwiseSqrt();
    Logger::Info(fmt::format("{}: converged={} ({} iters); position error {:.3f} m (1-sigma "
                             "{:.3f} m), velocity error {:.6f} m/s (1-sigma {:.6f} m/s)",
                             name_, converged_, num_iterations_,
                             (x0_est_.head(3) - x0_true_.head(3)).norm(), sigma.head(3).norm(),
                             (x0_est_.tail(3) - x0_true_.tail(3)).norm(), sigma.tail(3).norm()),
                 "GroundStationManagerApp");
    DataLogger::Log(fmt::format("{}/x0_true", name_), x0_true_);
    DataLogger::Log(fmt::format("{}/x0_estimated", name_), x0_est_);
    DataLogger::Log(fmt::format("{}/x0_sigma", name_), sigma);

    if (!run_srif_) return;

    // ---- Square-Root Information Filter + Dyer--McReynolds smoother ----
    const int n_state = 6;
    std::vector<std::vector<int>> obs_by_epoch(n_epochs);
    for (int k = 0; k < n_meas; ++k) obs_by_epoch[meas_[k].epoch_index].push_back(k);

    SRIF srif;
    srif.SetName(name_ + "/SRIF");
    srif.SetDynamicsFunction(
        [&](const State& x, Real t0, Real tf, const State* /*u*/, MatXd* F) -> State {
          MatXd stm;
          State xf = srif_dynamics_->Propagate(x, t0, tf, nullptr, &stm);
          if (F) *F = stm;
          return xf;
        });
    const bool use_pn = srif_use_process_noise_;
    const double psd = srif_accel_psd_;
    srif.SetProcessNoiseFunction([&](const State& /*x*/, Real t0, Real tf) -> MatXd {
      if (!use_pn || psd <= 0.0) return MatXd::Zero(n_state, n_state);
      return CwnaProcessNoise((tf - t0).val(), psd);
    });

    srif.SetState(Cart6(x0_est_.cast<Real>(), srif_dynamics_->GetFrame()));
    MatXd P0 = MatXd::Zero(n_state, n_state);
    for (int j = 0; j < 3; ++j) P0(j, j) = initial_position_sigma_m_ * initial_position_sigma_m_;
    for (int j = 3; j < 6; ++j)
      P0(j, j) = initial_velocity_sigma_mps_ * initial_velocity_sigma_mps_;
    srif.SetCovariance(P0);
    srif.SetTime(t_grid_(0));
    srif.InitializeLogger(n_epochs);

    srif_filtered_state_.resize(n_epochs, 6);
    srif_filtered_cov_.resize(n_epochs, 36);
    srif_smoothed_state_.resize(n_epochs, 6);
    srif_smoothed_cov_.resize(n_epochs, 36);

    for (int i = 0; i < n_epochs; ++i) {
      if (i > 0) srif.Predict(t_grid_(i));
      const std::vector<int>& obs_i = obs_by_epoch[i];

      srif.SetMeasurementFunction([&, i](const State& x, MatXd* H, MatXd* R) -> VecXd {
        Vec6d xm = x.cast<double>();
        // Row count for this epoch across its visible stations.
        int M = 0;
        for (int k : obs_i) M += (meas_[k].has_range ? 1 : 0) + (meas_[k].has_range_rate ? 1 : 0);
        VecXd z(M);
        MatXd Hm = MatXd::Zero(M, n_state);
        MatXd Rm = MatXd::Zero(M, M);
        int row = 0;
        for (int k : obs_i) {
          const StationMeasurement& m = meas_[k];
          MatXd h_obs;
          VecXd y = Observe(m, xm, &h_obs);
          for (int r = 0; r < y.size(); ++r) {
            z(row) = y(r);
            Hm.row(row) = h_obs.row(r);
            double s_range = RangeSigmaEffective(m);
            Rm(row, row) = (r == 0 && m.has_range) ? s_range * s_range
                                                   : m.range_rate_sigma * m.range_rate_sigma;
            row++;
          }
        }
        if (H) *H = Hm;
        if (R) *R = Rm;
        return z;
      });

      int M = 0;
      for (int k : obs_i) M += (meas_[k].has_range ? 1 : 0) + (meas_[k].has_range_rate ? 1 : 0);
      VecX z_true(M);
      int row = 0;
      for (int k : obs_i) {
        const StationMeasurement& m = meas_[k];
        if (m.has_range) z_true(row++) = m.range;
        if (m.has_range_rate) z_true(row++) = m.range_rate;
      }
      srif.Update(z_true);
      srif.LogFilterEstimate(i);
      srif_filtered_state_.row(i) = srif.GetState().cast<double>().transpose();
      Flatten6x6(srif.GetCovariance(), srif_filtered_cov_, i);
    }

    srif.InitializeSmootherState();
    srif_smoothed_state_.row(n_epochs - 1) = srif.GetSmoothedState(n_epochs - 1).transpose();
    Flatten6x6(srif.GetSmoothedCovariance(n_epochs - 1), srif_smoothed_cov_, n_epochs - 1);
    for (int i = n_epochs - 2; i >= 0; --i) {
      srif.UpdateSmoother(i);
      srif_smoothed_state_.row(i) = srif.GetSmoothedState(i).transpose();
      Flatten6x6(srif.GetSmoothedCovariance(i), srif_smoothed_cov_, i);
    }
  }

  void GroundStationManagerApp::Log(Real /*t*/) {
    if (!solved_) return;
    DataLogger::Log(fmt::format("{}/converged", name_), converged_ ? 1.0 : 0.0);
    DataLogger::Log(fmt::format("{}/num_measurements", name_), static_cast<double>(meas_.size()));
  }

  REGISTER_FACTORY_CLASS(Application, GroundStationManagerApp)

}  // namespace lupnt
