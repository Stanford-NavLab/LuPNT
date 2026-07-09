#include "lupnt/simulations/ground_station_odts/ground_station_odts_simulation.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <utility>

#include "lupnt/lupnt.h"
#include "lupnt/numerics/filters/srif.h"

namespace lupnt {

  std::vector<GroundStationOdtsStationConfig> DefaultDsnStations() {
    return {
        {"DSS14 (Goldstone)", 35.426456, 243.110461, 1001.39},
        {"DSS43 (Canberra)", -35.402424, 148.981267, 688.867},
        {"DSS63 (Madrid)", 40.431210, 355.751991, 864.816},
    };
  }

  namespace {
    Ptr<NBodyDynamics> BuildDynamics(const GroundStationOdtsConfig& cfg) {
      auto dynamics = MakePtr<NBodyDynamics>();
      dynamics->SetIntegrator(IntegratorType::RKF45);
      // abstol is a floor in meters; 1e-12 m is below double precision's representable floor at
      // ~1e6 m position magnitudes (~1e-10 m) and can make the adaptive step-size control
      // pathological. Match the tolerances validated for this ELFO scenario in Example 1/5/7.
      dynamics->SetIntegratorParams(IntegratorParams(1000, 1.0e-6, 1.0e-10));
      dynamics->AddBody(Body::Moon(cfg.moon_gravity_degree, cfg.moon_gravity_order));
      if (cfg.include_earth) dynamics->AddBody(Body::Earth());
      if (cfg.include_sun) dynamics->AddBody(Body::Sun());
      dynamics->SetFrame(Frame::MOON_CI);
      dynamics->SetTimeStep(cfg.integration_step_s);
      // Enable autodiff so Propagate(..., stm) can return the analytic
      // state-transition matrix used by the batch filter's analytic Jacobian and
      // by the covariance propagation. Value-only propagation is unaffected.
      dynamics->SetAutodiff(true);
      dynamics->SetUseRelativity(cfg.use_relativity);
      if (cfg.use_srp) dynamics->SetSrpCoefficient(cfg.srp_cr, cfg.srp_area_m2, cfg.srp_mass_kg);
      return dynamics;
    }

    int NumEpochs(const GroundStationOdtsConfig& cfg) {
      return static_cast<int>(std::floor(cfg.duration_s / cfg.obs_interval_s)) + 1;
    }

    int NumMeasurementRows(const GroundStationOdtsConfig& cfg) {
      return (cfg.use_range ? 1 : 0) + (cfg.use_range_rate ? 1 : 0);
    }
  }  // namespace

  GroundStationOdtsSimulation::GroundStationOdtsSimulation(GroundStationOdtsConfig config)
      : config_(std::move(config)) {}

  void GroundStationOdtsSimulation::Setup() {
    LUPNT_CHECK(!config_.ground_stations.empty(),
                "GroundStationODTS requires at least one ground station", "GroundStationODTS");
    LUPNT_CHECK(config_.use_range || config_.use_range_rate,
                "GroundStationODTS requires at least one measurement type (range and/or "
                "range-rate)",
                "GroundStationODTS");

    t0_tdb_ = ConvertTime(GregorianToTime(config_.start_epoch_utc), Time::UTC, Time::TDB);

    const int n_epochs = NumEpochs(config_);
    t_tdb_grid_.resize(n_epochs);
    results_.t_tdb.resize(n_epochs);
    for (int i = 0; i < n_epochs; ++i) {
      Real t = t0_tdb_ + static_cast<double>(i) * config_.obs_interval_s;
      t_tdb_grid_(i) = t;
      results_.t_tdb(i) = t.val();
    }

    BodyData earth = GetBodyData(BodyId::EARTH);
    results_.station_names.clear();
    station_r_ecef_.clear();
    for (const auto& station : config_.ground_stations) {
      results_.station_names.push_back(station.name);

      LatLonAlt lla(Vec3(station.latitude_deg, station.longitude_deg, station.altitude_m),
                    earth.fixed_frame);
      Vec3 r_ecef = LatLonAltToCart(lla, earth.R, earth.flattening);
      station_r_ecef_.push_back(r_ecef);
    }

    setup_complete_ = true;
    precompute_complete_ = false;
  }

  void GroundStationOdtsSimulation::Precompute() {
    if (!setup_complete_) Setup();

    const int n_epochs = static_cast<int>(t_tdb_grid_.size());
    const int n_stations = static_cast<int>(config_.ground_stations.size());

    // Truth orbit: ELFO classical elements in the Moon-centered Orbital-Plane frame,
    // converted to the Moon-Centered Inertial integration frame.
    Vec6 coe_op(config_.orbit_a_m, config_.orbit_ecc, config_.orbit_inc_rad, config_.orbit_raan_rad,
                config_.orbit_argp_rad, config_.orbit_mean_anomaly_rad);
    Vec6 rv0_op = ClassicalToCart(coe_op, GM_MOON);
    Vec6 x0_true = ConvertFrame(t0_tdb_, rv0_op, Frame::MOON_OP, Frame::MOON_CI);
    results_.x0_true = x0_true.cast<double>();

    auto dynamics = BuildDynamics(config_);
    MatX6 truth_mci(n_epochs, 6);
    {
      State x = Cart6(x0_true, Frame::MOON_CI);
      truth_mci.row(0) = x.transpose();
      for (int i = 1; i < n_epochs; ++i) {
        x = dynamics->Propagate(x, t_tdb_grid_(i - 1), t_tdb_grid_(i), nullptr);
        truth_mci.row(i) = x.transpose();
      }
    }
    results_.truth_state = truth_mci.cast<double>();

    // Visibility: topocentric elevation of the truth trajectory at each station.
    MatX6 truth_ecef = ConvertFrame(t_tdb_grid_, truth_mci, Frame::MOON_CI, Frame::ECEF);
    results_.elevation_deg.resize(n_epochs, n_stations);
    for (int s = 0; s < n_stations; ++s) {
      Cart3 r_gs(station_r_ecef_[s], Frame::ECEF);
      for (int i = 0; i < n_epochs; ++i) {
        Cart3 r_sat(Vec3(truth_ecef.row(i).head(3).transpose()), Frame::ECEF);
        State aer = CartToAzElRange(r_sat, r_gs);
        results_.elevation_deg(i, s) = (aer(1) * DEG).val();
      }
    }

    // Simulated measurements: range/range-rate + noise at every (epoch, station) pair
    // above the elevation mask.
    std::mt19937 noise_rng(static_cast<unsigned int>(config_.seed));
    std::vector<int> obs_epoch, obs_station;
    std::vector<double> obs_rho_true, obs_rhodot_true, obs_rho, obs_rhodot;

    for (int i = 0; i < n_epochs; ++i) {
      Vec3 r_sat = truth_ecef.row(i).head(3).transpose();
      Vec3 v_sat = truth_ecef.row(i).tail(3).transpose();
      for (int s = 0; s < n_stations; ++s) {
        if (results_.elevation_deg(i, s) <= config_.elevation_mask_deg) continue;

        Vec3 dr = r_sat - station_r_ecef_[s];
        Real rho = dr.norm();
        Real rhodot = dr.dot(v_sat) / rho;

        obs_epoch.push_back(i);
        obs_station.push_back(s);
        obs_rho_true.push_back(rho.val());
        obs_rhodot_true.push_back(rhodot.val());
        obs_rho.push_back(config_.use_range
                              ? rho.val()
                                    + SampleNormal(0.0, config_.range_sigma_m, &noise_rng).val()
                              : std::numeric_limits<double>::quiet_NaN());
        obs_rhodot.push_back(
            config_.use_range_rate
                ? rhodot.val() + SampleNormal(0.0, config_.range_rate_sigma_mps, &noise_rng).val()
                : std::numeric_limits<double>::quiet_NaN());
      }
    }

    const int n_obs = static_cast<int>(obs_epoch.size());
    LUPNT_CHECK(n_obs > 0,
                "No station visibility above the elevation mask over the simulated arc -- widen "
                "duration_s or relax elevation_mask_deg",
                "GroundStationODTS");

    results_.obs_epoch_index.resize(n_obs);
    results_.obs_station_index.resize(n_obs);
    results_.obs_range_true_m.resize(n_obs);
    results_.obs_range_rate_true_mps.resize(n_obs);
    results_.obs_range_m.resize(n_obs);
    results_.obs_range_rate_mps.resize(n_obs);
    for (int k = 0; k < n_obs; ++k) {
      results_.obs_epoch_index(k) = obs_epoch[k];
      results_.obs_station_index(k) = obs_station[k];
      results_.obs_range_true_m(k) = obs_rho_true[k];
      results_.obs_range_rate_true_mps(k) = obs_rhodot_true[k];
      results_.obs_range_m(k) = obs_rho[k];
      results_.obs_range_rate_mps(k) = obs_rhodot[k];
    }

    precompute_complete_ = true;
  }

  void GroundStationOdtsSimulation::Run() {
    if (!precompute_complete_) Precompute();

    const int n_epochs = static_cast<int>(t_tdb_grid_.size());
    const int n_stations = static_cast<int>(config_.ground_stations.size());
    const int n_obs = static_cast<int>(results_.obs_epoch_index.size());
    const int n_rows = NumMeasurementRows(config_);

    // A-priori state error injected into the batch filter's starting guess.
    std::mt19937 guess_rng(static_cast<unsigned int>(config_.seed) + 1);
    VecXd dx0(6);
    for (int j = 0; j < 3; ++j) {
      dx0(j) = SampleNormal(0.0, config_.initial_position_sigma_m, &guess_rng).val();
    }
    for (int j = 3; j < 6; ++j) {
      dx0(j) = SampleNormal(0.0, config_.initial_velocity_sigma_mps, &guess_rng).val();
    }
    results_.x0_initial_guess = results_.x0_true + dx0;

    VecXd initial_covariance_diag(6);
    initial_covariance_diag.head(3).setConstant(config_.initial_position_sigma_m
                                                * config_.initial_position_sigma_m);
    initial_covariance_diag.tail(3).setConstant(config_.initial_velocity_sigma_mps
                                                * config_.initial_velocity_sigma_mps);

    // Measurement vectors/weights, in the same (epoch, station) order as Precompute()
    // populated results_.obs_*.
    const double w_range = 1.0 / (config_.range_sigma_m * config_.range_sigma_m);
    const double w_range_rate = 1.0 / (config_.range_rate_sigma_mps * config_.range_rate_sigma_mps);
    std::vector<VecXd> measurements;
    std::vector<VecXd> measurement_weights;
    measurements.reserve(n_obs);
    measurement_weights.reserve(n_obs);
    for (int k = 0; k < n_obs; ++k) {
      VecXd y(n_rows);
      VecXd w(n_rows);
      int row = 0;
      if (config_.use_range) {
        y(row) = results_.obs_range_m(k);
        w(row) = w_range;
        row++;
      }
      if (config_.use_range_rate) {
        y(row) = results_.obs_range_rate_mps(k);
        w(row) = w_range_rate;
        row++;
      }
      measurements.push_back(y);
      measurement_weights.push_back(w);
    }

    auto dynamics_filter = BuildDynamics(config_);
    const std::vector<Vec3>& station_r_ecef = station_r_ecef_;
    const VecXi& obs_epoch_index = results_.obs_epoch_index;
    const VecXi& obs_station_index = results_.obs_station_index;
    const VecX& t_tdb_grid = t_tdb_grid_;
    const GroundStationOdtsConfig& cfg = config_;

    // Station states in the MOON_CI (inertial) frame at every epoch: the antenna is
    // fixed in ECEF (zero ECEF velocity), so its MOON_CI velocity is purely the
    // frame-rotation term. Needed for the analytic range/range-rate partials below,
    // which treat the station state as a known, x0-independent quantity.
    std::vector<MatXd> station_state_mci(n_stations);
    for (int s = 0; s < n_stations; ++s) {
      MatX6 gs_ecef = MatX6::Zero(n_epochs, 6);
      for (int i = 0; i < n_epochs; ++i) gs_ecef.row(i).head(3) = station_r_ecef_[s].transpose();
      station_state_mci[s]
          = ConvertFrame(t_tdb_grid_, gs_ecef, Frame::ECEF, Frame::MOON_CI).cast<double>();
    }

    // Propagate an epoch state x0 over the full observation grid, returning the
    // MOON_CI trajectory and the cumulative state-transition matrices
    // Phi(t_i, t0) = d x(t_i) / d x0. Each per-segment STM d x(t_i)/d x(t_{i-1}) is the
    // autodiff (analytic) Jacobian from NBodyDynamics::Propagate(..., stm); the segment
    // STMs are chained to give Phi(t_i, t0). No finite differencing is used.
    auto PropagateGridStm = [&](const VecXd& x0, MatX6& grid_mci, std::vector<Mat6d>& stm_cum) {
      grid_mci.resize(n_epochs, 6);
      stm_cum.assign(n_epochs, Mat6d::Identity());
      State x = Cart6(x0.cast<Real>(), Frame::MOON_CI);
      grid_mci.row(0) = x.transpose();
      Mat6d phi = Mat6d::Identity();
      for (int i = 1; i < n_epochs; ++i) {
        MatXd stm_seg;
        x = dynamics_filter->Propagate(x, t_tdb_grid(i - 1), t_tdb_grid(i), nullptr, &stm_seg);
        phi = stm_seg * phi;
        stm_cum[i] = phi;
        grid_mci.row(i) = x.transpose();
      }
    };

    // Batch design matrix d(range, range_rate)/d(x0). The predicted measurements are
    // evaluated in ECEF via the identical rho/rho_dot geometry Precompute() used for the
    // "truth" measurements (dr = r_sat_ecef - station_r_ecef, rho_dot = dr.dot(v_sat_ecef)/rho,
    // station velocity zero in ECEF), so the filter's predictions share the observed
    // measurements' geometric pathway. The Jacobian, however, is formed analytically: the
    // closed-form range/range-rate partials w.r.t. the inertial (MOON_CI) satellite state --
    // range and range-rate are frame-invariant, so they are evaluated with the station's
    // MOON_CI state -- chained with the autodiff STM Phi(t_i, t0). The cache is rebuilt
    // whenever meas_idx == 0, i.e. once per RunBatchFilter() pass -- see batch_filter.cc.
    struct AnalyticCache {
      MatX6 grid_mci;              // [n_epochs x 6]
      MatXd grid_ecef;             // [n_epochs x 6]
      std::vector<Mat6d> stm_cum;  // [n_epochs] of Phi(t_i, t0)
    };
    AnalyticCache acache;

    MeasurementModelFunction measurement_model_analytic
        = [&](const VecXd& state_estimate, int meas_idx) -> std::pair<VecXd, MatXd> {
      if (meas_idx == 0) {
        PropagateGridStm(state_estimate, acache.grid_mci, acache.stm_cum);
        acache.grid_ecef
            = ConvertFrame(t_tdb_grid, acache.grid_mci, Frame::MOON_CI, Frame::ECEF).cast<double>();
      }

      const int i = obs_epoch_index(meas_idx);
      const int s = obs_station_index(meas_idx);

      // Predicted measurement in ECEF (station fixed, zero ECEF velocity).
      Vec3d r_gs_e = station_r_ecef[s].cast<double>();
      Vec3d r_sat_e = acache.grid_ecef.row(i).head(3).transpose();
      Vec3d v_sat_e = acache.grid_ecef.row(i).tail(3).transpose();
      Vec3d dr_e = r_sat_e - r_gs_e;
      double rho = dr_e.norm();
      double rho_dot = dr_e.dot(v_sat_e) / rho;

      // Analytic observation partials w.r.t. the MOON_CI satellite state. range and
      // range-rate are scalar (frame-invariant), evaluated here with the station's MOON_CI
      // state so the partials chain directly with the MOON_CI STM.
      Vec6d xs_m = acache.grid_mci.row(i).cast<double>().transpose();
      Vec6d xgs_m = station_state_mci[s].row(i).transpose();
      Vec3d dr_m = xs_m.head(3) - xgs_m.head(3);
      Vec3d dv_m = xs_m.tail(3) - xgs_m.tail(3);
      double rho_m = dr_m.norm();
      Vec3d u = dr_m / rho_m;  // d(rho)/d(r)
      double rho_dot_m = dr_m.dot(dv_m) / rho_m;

      MatXd h_obs = MatXd::Zero(n_rows, 6);
      VecXd y(n_rows);
      int row = 0;
      if (cfg.use_range) {
        h_obs.block(row, 0, 1, 3) = u.transpose();
        y(row) = rho;
        row++;
      }
      if (cfg.use_range_rate) {
        h_obs.block(row, 0, 1, 3) = ((dv_m - rho_dot_m * u) / rho_m).transpose();
        h_obs.block(row, 3, 1, 3) = u.transpose();
        y(row) = rho_dot;
        row++;
      }
      MatXd h_local = h_obs * acache.stm_cum[i];
      return {y, h_local};
    };

    // Legacy finite-difference (numerical) Jacobian, retained for comparison/fallback
    // behind cfg.batch_use_analytic_jacobian. Propagates x0 and 6 perturbed copies over
    // the grid and differences the ECEF observation model.
    struct NumericCache {
      MatXd state_ecef_nominal;              // [n_epochs x 6]
      std::array<MatXd, 6> state_ecef_pert;  // [6] of [n_epochs x 6]
    };
    auto ncache = std::make_shared<NumericCache>();
    const double eps_pos = 1.0;     // [m]   finite-difference step, position
    const double eps_vel = 1.0e-3;  // [m/s] finite-difference step, velocity
    auto PropagateFullGridEcef = [&](const VecXd& x0) -> MatXd {
      MatX6 grid_mci(n_epochs, 6);
      State x = Cart6(x0.cast<Real>(), Frame::MOON_CI);
      grid_mci.row(0) = x.transpose();
      for (int i = 1; i < n_epochs; ++i) {
        x = dynamics_filter->Propagate(x, t_tdb_grid(i - 1), t_tdb_grid(i), nullptr);
        grid_mci.row(i) = x.transpose();
      }
      return ConvertFrame(t_tdb_grid, grid_mci, Frame::MOON_CI, Frame::ECEF).cast<double>();
    };
    MeasurementModelFunction measurement_model_numeric
        = [&](const VecXd& state_estimate, int meas_idx) -> std::pair<VecXd, MatXd> {
      if (meas_idx == 0) {
        ncache->state_ecef_nominal = PropagateFullGridEcef(state_estimate);
        for (int d = 0; d < 6; ++d) {
          VecXd x0_pert = state_estimate;
          x0_pert(d) += (d < 3 ? eps_pos : eps_vel);
          ncache->state_ecef_pert[d] = PropagateFullGridEcef(x0_pert);
        }
      }
      const int epoch_idx = obs_epoch_index(meas_idx);
      const int station_idx = obs_station_index(meas_idx);
      Vec3d r_gs_d = station_r_ecef[station_idx].cast<double>();
      auto ComputeObs = [&](const MatXd& grid_ecef) -> VecXd {
        Vec3d r_sat = grid_ecef.row(epoch_idx).head(3).transpose();
        Vec3d v_sat = grid_ecef.row(epoch_idx).tail(3).transpose();
        Vec3d dr = r_sat - r_gs_d;
        double rho = dr.norm();
        double rho_dot = dr.dot(v_sat) / rho;
        VecXd y_obs(n_rows);
        int row = 0;
        if (cfg.use_range) {
          y_obs(row) = rho;
          row++;
        }
        if (cfg.use_range_rate) {
          y_obs(row) = rho_dot;
          row++;
        }
        return y_obs;
      };
      VecXd y = ComputeObs(ncache->state_ecef_nominal);
      MatXd h_local(n_rows, 6);
      for (int d = 0; d < 6; ++d) {
        VecXd y_pert = ComputeObs(ncache->state_ecef_pert[d]);
        h_local.col(d) = (y_pert - y) / (d < 3 ? eps_pos : eps_vel);
      }
      return {y, h_local};
    };

    MeasurementModelFunction measurement_model = config_.batch_use_analytic_jacobian
                                                     ? measurement_model_analytic
                                                     : measurement_model_numeric;

    BatchFilterConfig bf_config;
    bf_config.use_weights = config_.batch_use_weights;
    bf_config.use_initialization = config_.batch_use_initialization;
    bf_config.convergence_tol = config_.batch_convergence_tol;
    bf_config.max_iterations = config_.batch_max_iterations;

    BatchFilterResults bf_results
        = RunBatchFilter(results_.x0_initial_guess, initial_covariance_diag, measurements,
                         measurement_weights, measurement_model, bf_config, results_.x0_true);

    results_.x0_estimated = bf_results.state_estimate;
    results_.covariance = bf_results.state_covariance;
    results_.converged = bf_results.converged;
    results_.num_iterations = bf_results.iterations;

    // Propagate the truth/estimated trajectories once more (value-only, no STM) to fill in the
    // per-iteration diagnostics and the full estimated-trajectory time series. Returns the
    // MOON_CI trajectory; callers convert to ECEF themselves when comparing against measurements
    // (see the rho/rho_dot comment on `measurement_model` above).
    auto PropagateGrid = [&](const VecXd& x0) -> MatX6 {
      auto dyn = BuildDynamics(config_);
      MatX6 grid(n_epochs, 6);
      State x = Cart6(x0.cast<Real>(), Frame::MOON_CI);
      grid.row(0) = x.transpose();
      for (int i = 1; i < n_epochs; ++i) {
        x = dyn->Propagate(x, t_tdb_grid_(i - 1), t_tdb_grid_(i), nullptr);
        grid.row(i) = x.transpose();
      }
      return grid;
    };

    const int n_iter = static_cast<int>(bf_results.iteration_history.size());
    results_.iteration_state_estimate.resize(n_iter, 6);
    results_.iteration_correction_norm.resize(n_iter);
    results_.iteration_weighted_rms.resize(n_iter);
    results_.iteration_pos_error_m.resize(n_iter);
    results_.iteration_vel_error_mps.resize(n_iter);
    results_.iteration_rms_range_m.resize(n_iter);
    results_.iteration_rms_range_rate_mps.resize(n_iter);

    for (int k = 0; k < n_iter; ++k) {
      const auto& info = bf_results.iteration_history[k];
      results_.iteration_state_estimate.row(k) = info.state_estimate.transpose();
      results_.iteration_correction_norm(k) = info.correction_norm;
      results_.iteration_weighted_rms(k) = info.weighted_rms;
      results_.iteration_pos_error_m(k)
          = (info.state_estimate.head(3) - results_.x0_true.head(3)).norm();
      results_.iteration_vel_error_mps(k)
          = (info.state_estimate.tail(3) - results_.x0_true.tail(3)).norm();

      MatXd grid_ecef = ConvertFrame(t_tdb_grid_, PropagateGrid(info.state_estimate),
                                     Frame::MOON_CI, Frame::ECEF)
                            .cast<double>();
      double sse_range = 0.0, sse_rate = 0.0;
      for (int k_obs = 0; k_obs < n_obs; ++k_obs) {
        int i = results_.obs_epoch_index(k_obs);
        int s = results_.obs_station_index(k_obs);
        Vec3d r_gs_d = station_r_ecef_[s].cast<double>();
        Vec3d r_rel = grid_ecef.row(i).head(3).transpose() - r_gs_d;
        Vec3d v_sat = grid_ecef.row(i).tail(3).transpose();
        double rho = r_rel.norm();
        double rho_dot = r_rel.dot(v_sat) / rho;
        if (config_.use_range) sse_range += std::pow(results_.obs_range_m(k_obs) - rho, 2);
        if (config_.use_range_rate) {
          sse_rate += std::pow(results_.obs_range_rate_mps(k_obs) - rho_dot, 2);
        }
      }
      results_.iteration_rms_range_m(k) = config_.use_range
                                              ? std::sqrt(sse_range / n_obs)
                                              : std::numeric_limits<double>::quiet_NaN();
      results_.iteration_rms_range_rate_mps(k) = config_.use_range_rate
                                                     ? std::sqrt(sse_rate / n_obs)
                                                     : std::numeric_limits<double>::quiet_NaN();
    }

    // Final estimated trajectory and its formal covariance propagated over the full grid.
    // Re-propagate x0_estimated with the autodiff STM Phi(t_i, t0) and map the epoch
    // covariance forward: P(t_i) = Phi(t_i, t0) * covariance * Phi(t_i, t0)^T.
    {
      MatX6 est_mci;
      std::vector<Mat6d> stm_cum;
      PropagateGridStm(results_.x0_estimated, est_mci, stm_cum);
      results_.estimated_state = est_mci.cast<double>();

      Mat6d P0 = results_.covariance;
      results_.estimated_covariance.resize(n_epochs, 36);
      for (int i = 0; i < n_epochs; ++i) {
        Mat6d Pi = stm_cum[i] * P0 * stm_cum[i].transpose();
        for (int r = 0; r < 6; ++r) {
          for (int c = 0; c < 6; ++c) results_.estimated_covariance(i, r * 6 + c) = Pi(r, c);
        }
      }
    }

    // Square-Root Information Filter (SRIF) + smoother. The `SRIF` filter (an EKF that
    // factors the information matrix, see lupnt/numerics/filters/srif.h) is driven over
    // the observation grid: it starts from the batch epoch solution with the same
    // diffuse a-priori as the batch, re-linearizes and processes each epoch's tracking
    // once forward (numerically-robust square-root form) with an injected
    // white-noise-acceleration process noise, then smooths backward with the inherited
    // RTS pass. The per-epoch smoothed covariance reflects all measurements plus the
    // process noise -- a more realistic (less optimistic) uncertainty than mapping the
    // batch epoch covariance forward with Phi P Phi^T.
    if (config_.run_srif) {
      const int n_state = 6;
      const int meas_rows = NumMeasurementRows(config_);
      const double range_var = config_.range_sigma_m * config_.range_sigma_m;
      const double rate_var = config_.range_rate_sigma_mps * config_.range_rate_sigma_mps;
      const bool srif_use_pn = config_.srif_use_process_noise;
      const double srif_psd = config_.srif_accel_psd;

      // Group observation indices by epoch (obs are stored in epoch order already).
      std::vector<std::vector<int>> obs_by_epoch(n_epochs);
      for (int k = 0; k < n_obs; ++k) obs_by_epoch[results_.obs_epoch_index(k)].push_back(k);

      SRIF srif;
      srif.SetName("SRIF");

      // Dynamics: propagate the MOON_CI Cartesian state and return the STM.
      srif.SetDynamicsFunction(
          [&](const State& x, Real t0, Real tf, const State* /*u*/, MatXd* F) -> State {
            MatXd stm;
            State xf = dynamics_filter->Propagate(x, t0, tf, nullptr, &stm);
            if (F) *F = stm;
            return xf;
          });

      // Continuous white-noise-acceleration process noise (models unmodeled dynamics).
      srif.SetProcessNoiseFunction([&](const State& /*x*/, Real t0, Real tf) -> MatXd {
        if (!srif_use_pn || srif_psd <= 0.0) return MatXd::Zero(n_state, n_state);
        return CwnaProcessNoise((tf - t0).val(), srif_psd);
      });

      // Same diffuse a-priori as the batch, started at the batch epoch solution.
      srif.SetState(Cart6(results_.x0_estimated.cast<Real>(), Frame::MOON_CI));
      MatXd P0 = MatXd::Zero(n_state, n_state);
      const double pos_var = config_.initial_position_sigma_m * config_.initial_position_sigma_m;
      const double vel_var
          = config_.initial_velocity_sigma_mps * config_.initial_velocity_sigma_mps;
      for (int j = 0; j < 3; ++j) P0(j, j) = pos_var;
      for (int j = 3; j < 6; ++j) P0(j, j) = vel_var;
      srif.SetCovariance(P0);
      srif.SetTime(t_tdb_grid_(0));
      srif.InitializeLogger(n_epochs);

      results_.srif_filtered_state.resize(n_epochs, 6);
      results_.srif_filtered_covariance.resize(n_epochs, 36);
      results_.srif_smoothed_state.resize(n_epochs, 6);
      results_.srif_smoothed_covariance.resize(n_epochs, 36);

      auto flatten = [](const MatXd& P, MatXd& dst, int i) {
        for (int r = 0; r < 6; ++r)
          for (int c = 0; c < 6; ++c) dst(i, r * 6 + c) = P(r, c);
      };

      // Forward pass: predict to each epoch, then update with that epoch's tracking.
      for (int i = 0; i < n_epochs; ++i) {
        if (i > 0) srif.Predict(t_tdb_grid_(i));

        const std::vector<int>& obs_i = obs_by_epoch[i];
        // Measurement model for epoch i: stacked range/range-rate for the visible
        // stations, predicted and linearized in MOON_CI. Range and range-rate are
        // frame-invariant, so this matches the ECEF geometry used to simulate the
        // observations, while the station's MOON_CI state carries the frame-rotation
        // velocity the range-rate partials need.
        srif.SetMeasurementFunction([&, i](const State& x, MatXd* H, MatXd* R) -> VecXd {
          Vec6d xm = x.cast<double>();
          const int M = meas_rows * static_cast<int>(obs_i.size());
          VecXd z(M);
          MatXd Hm = MatXd::Zero(M, n_state);
          MatXd Rm = MatXd::Zero(M, M);
          int row = 0;
          for (int k : obs_i) {
            const int s = results_.obs_station_index(k);
            Vec6d xgs = station_state_mci[s].row(i).transpose();
            Vec3d dr = xm.head(3) - xgs.head(3);
            Vec3d dv = xm.tail(3) - xgs.tail(3);
            double rho = dr.norm();
            Vec3d uhat = dr / rho;
            double rho_dot = dr.dot(dv) / rho;
            if (config_.use_range) {
              z(row) = rho;
              Hm.block(row, 0, 1, 3) = uhat.transpose();
              Rm(row, row) = range_var;
              row++;
            }
            if (config_.use_range_rate) {
              z(row) = rho_dot;
              Hm.block(row, 0, 1, 3) = ((dv - rho_dot * uhat) / rho).transpose();
              Hm.block(row, 3, 1, 3) = uhat.transpose();
              Rm(row, row) = rate_var;
              row++;
            }
          }
          if (H) *H = Hm;
          if (R) *R = Rm;
          return z;
        });

        // Observed measurement vector for epoch i, in the same station/row order.
        const int M = meas_rows * static_cast<int>(obs_i.size());
        VecX z_true(M);
        int row = 0;
        for (int k : obs_i) {
          if (config_.use_range) z_true(row++) = results_.obs_range_m(k);
          if (config_.use_range_rate) z_true(row++) = results_.obs_range_rate_mps(k);
        }
        srif.Update(z_true);
        srif.LogFilterEstimate(i);

        results_.srif_filtered_state.row(i) = srif.GetState().cast<double>().transpose();
        flatten(srif.GetCovariance(), results_.srif_filtered_covariance, i);
      }

      // Backward pass: fixed-interval RTS smoother (inherited from EKF).
      srif.InitializeSmootherState();
      results_.srif_smoothed_state.row(n_epochs - 1)
          = srif.GetSmoothedState(n_epochs - 1).transpose();
      flatten(srif.GetSmoothedCovariance(n_epochs - 1), results_.srif_smoothed_covariance,
              n_epochs - 1);
      for (int i = n_epochs - 2; i >= 0; --i) {
        srif.UpdateSmoother(i);
        results_.srif_smoothed_state.row(i) = srif.GetSmoothedState(i).transpose();
        flatten(srif.GetSmoothedCovariance(i), results_.srif_smoothed_covariance, i);
      }
    }
  }

}  // namespace lupnt
