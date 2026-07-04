#include "lupnt/simulations/GroundStationODTS/ground_station_odts_simulation.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <utility>

#include "lupnt/lupnt.h"

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
    Vec6 coe_op(config_.orbit_a_m, config_.orbit_ecc, config_.orbit_inc_rad,
               config_.orbit_raan_rad, config_.orbit_argp_rad, config_.orbit_mean_anomaly_rad);
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
                              ? rho.val() + SampleNormal(0.0, config_.range_sigma_m, &noise_rng).val()
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

    // Batch design matrix: propagate x0 (and 6 finite-difference-perturbed copies of x0)
    // sequentially over the observation grid with the plain (no-STM) NBodyDynamics::Propagate,
    // convert the resulting MOON_CI trajectory to ECEF, and evaluate the analytic range/
    // range-rate observation model there -- exactly mirroring the rho/rho_dot formula
    // Precompute() uses to generate the "truth" measurements (dr = r_sat_ecef - station_r_ecef,
    // rho_dot = dr.dot(v_sat_ecef)/rho, station velocity zero in ECEF) -- so the filter's
    // predicted measurements are computed via the identical geometric pathway as the observed
    // ones. Differentiating that model numerically gives d(range, range_rate)/d(x0). The cache
    // is rebuilt whenever meas_idx == 0, i.e. once per RunBatchFilter() pass over
    // `measurements` -- see batch_filter.cc.
    struct FilterCache {
      MatXd state_ecef_nominal;                 // [n_epochs x 6]
      std::array<MatXd, 6> state_ecef_pert;      // [6] of [n_epochs x 6]
    };
    auto cache = std::make_shared<FilterCache>();

    auto dynamics_filter = BuildDynamics(config_);
    const std::vector<Vec3>& station_r_ecef = station_r_ecef_;
    const VecXi& obs_epoch_index = results_.obs_epoch_index;
    const VecXi& obs_station_index = results_.obs_station_index;
    const VecX& t_tdb_grid = t_tdb_grid_;
    const GroundStationOdtsConfig& cfg = config_;

    const double eps_pos = 1.0;      // [m]   finite-difference step, position
    const double eps_vel = 1.0e-3;   // [m/s] finite-difference step, velocity

    auto PropagateFullGridEcef = [dynamics_filter, t_tdb_grid, n_epochs](const VecXd& x0) -> MatXd {
      MatX6 grid_mci(n_epochs, 6);
      State x = Cart6(x0.cast<Real>(), Frame::MOON_CI);
      grid_mci.row(0) = x.transpose();
      for (int i = 1; i < n_epochs; ++i) {
        x = dynamics_filter->Propagate(x, t_tdb_grid(i - 1), t_tdb_grid(i), nullptr);
        grid_mci.row(i) = x.transpose();
      }
      MatXd grid_ecef = ConvertFrame(t_tdb_grid, grid_mci, Frame::MOON_CI, Frame::ECEF).cast<double>();
      return grid_ecef;
    };

    MeasurementModelFunction measurement_model = [=](const VecXd& state_estimate,
                                                      int meas_idx) -> std::pair<VecXd, MatXd> {
      if (meas_idx == 0) {
        cache->state_ecef_nominal = PropagateFullGridEcef(state_estimate);
        for (int d = 0; d < 6; ++d) {
          VecXd x0_pert = state_estimate;
          x0_pert(d) += (d < 3 ? eps_pos : eps_vel);
          cache->state_ecef_pert[d] = PropagateFullGridEcef(x0_pert);
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

      VecXd y = ComputeObs(cache->state_ecef_nominal);
      MatXd h_local(n_rows, 6);
      for (int d = 0; d < 6; ++d) {
        VecXd y_pert = ComputeObs(cache->state_ecef_pert[d]);
        h_local.col(d) = (y_pert - y) / (d < 3 ? eps_pos : eps_vel);
      }
      return {y, h_local};
    };

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
      results_.iteration_pos_error_m(k) = (info.state_estimate.head(3) - results_.x0_true.head(3)).norm();
      results_.iteration_vel_error_mps(k)
          = (info.state_estimate.tail(3) - results_.x0_true.tail(3)).norm();

      MatXd grid_ecef
          = ConvertFrame(t_tdb_grid_, PropagateGrid(info.state_estimate), Frame::MOON_CI, Frame::ECEF)
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
      results_.iteration_rms_range_m(k)
          = config_.use_range ? std::sqrt(sse_range / n_obs) : std::numeric_limits<double>::quiet_NaN();
      results_.iteration_rms_range_rate_mps(k) = config_.use_range_rate
                                                     ? std::sqrt(sse_rate / n_obs)
                                                     : std::numeric_limits<double>::quiet_NaN();
    }

    results_.estimated_state = PropagateGrid(results_.x0_estimated).cast<double>();
  }

}  // namespace lupnt
