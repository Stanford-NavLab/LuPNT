#include "src/filter_setup.h"

#include <lupnt/lupnt.h>

#include <iomanip>

#include "src/delays.h"
#include "src/dynamics_models.h"
#include "src/measurement_storage.h"
#include "src/state_meas_manager.h"

namespace filtering_sim {
  using namespace lupnt;

  /***************************************************************************
   * Filter Setup Functions
   * **************************************************************************/

  bool AddMeasFlags(const YAML::Node& meas_config, const std::string& meas_type_str,
                    MeasType meas_type, std::map<FreqTypePair, MeasFlag>& meas_flags) {
    auto node = meas_config[meas_type_str];
    bool is_ionofree
        = (meas_type == MeasType::IONOFREE_CODE || meas_type == MeasType::IONOFREE_CARRIER
           || meas_type == MeasType::IONOFREE_TDCP);
    bool use_measurement = false;

    if (!node) {
      return use_measurement;
    }
    if (!is_ionofree) {
      std::vector<double> L1_config = node["L1"].as<std::vector<double>>();
      MeasFlag flag_L1;
      flag_L1.use_meas = static_cast<bool>(L1_config[0]);
      flag_L1.min_alt = static_cast<double>(L1_config[1]);
      flag_L1.max_alt = static_cast<double>(L1_config[2]);
      flag_L1.process_inv = static_cast<int>(L1_config[3]);
      use_measurement |= flag_L1.use_meas;
      std::vector<double> L5_config = node["L5"].as<std::vector<double>>();
      MeasFlag flag_L5;
      flag_L5.use_meas = static_cast<bool>(L5_config[0]);
      flag_L5.min_alt = static_cast<double>(L5_config[1]);
      flag_L5.max_alt = static_cast<double>(L5_config[2]);
      flag_L5.process_inv = static_cast<int>(L5_config[3]);
      use_measurement |= flag_L5.use_meas;
      meas_flags[std::make_pair(MeasFreq::L1, meas_type)] = flag_L1;
      meas_flags[std::make_pair(MeasFreq::L5, meas_type)] = flag_L5;
    } else {
      std::vector<double> ionofree_config = node["L1_L5"].as<std::vector<double>>();
      MeasFlag flag_ionofree;
      flag_ionofree.use_meas = static_cast<bool>(ionofree_config[0]);
      flag_ionofree.min_alt = static_cast<double>(ionofree_config[1]);
      flag_ionofree.max_alt = static_cast<double>(ionofree_config[2]);
      flag_ionofree.process_inv = static_cast<int>(ionofree_config[3]);
      use_measurement |= flag_ionofree.use_meas;
      meas_flags[std::make_pair(MeasFreq::L1_L5, meas_type)] = flag_ionofree;
    }

    return use_measurement;
  }

  FilterConfig SetupFilterConfig(const SimulationConfig& sim_config, const VecXd& tspan, int seed,
                                 bool verbose) {
    FilterConfig cfg;

    // Load filter parameters from sim_config
    auto filter_node = sim_config.config["filter"];

    // Process Noise
    // ---------------------------------------------------------------------------------
    auto process_noise_cfg = sim_config.config["process_noise"];
    int N_window = process_noise_cfg["N_window"].as<int>(10);
    int N_wait = process_noise_cfg["N_wait"].as<int>(0);
    double alpha = process_noise_cfg["alpha"].as<double>(0.1);
    MatXd Q_a = process_noise_cfg["sigma_a"].as<VecXd>().array().square().matrix().asDiagonal();
    Vec3d beta = process_noise_cfg["beta"].as<Vec3d>();
    double Q_sigma_min = process_noise_cfg["Q_sigma_min"].as<double>(1e-8);
    double Q_sigma_max = process_noise_cfg["Q_sigma_max"].as<double>(1e8);
    int use_rtn = process_noise_cfg["rtn"].as<bool>(false);

    auto alg_str = process_noise_cfg["algorithm"].as<std::string>();
    ProcessNoiseAlgorithm alg;
    if (alg_str == "SNC") {
      alg = ProcessNoiseAlgorithm::SNC;
    } else if (alg_str == "ASNC") {
      alg = ProcessNoiseAlgorithm::ASNC;
    } else if (alg_str == "ADMC") {
      alg = ProcessNoiseAlgorithm::ADMC;
    } else {
      throw std::runtime_error("Unknown process noise algorithm: " + alg_str);
    }

    ProcessNoise process_noise;
    process_noise.SetAlgorithm(alg);
    process_noise.SetNoiseLimits(Q_sigma_min, Q_sigma_max);
    process_noise.SetTimeStep(sim_config.dt);
    process_noise.SetWindowSize(N_window);
    process_noise.SetProcessNoise(Q_a);
    process_noise.SetBeta(beta);
    process_noise.SetAlpha(alpha);
    process_noise.SetWaitSteps(N_wait);

    cfg.Q_a = Q_a;
    cfg.process_noise = process_noise;
    cfg.process_noise_algorithm = alg;
    cfg.dmc = (alg == ProcessNoiseAlgorithm::ADMC);
    cfg.dmc_beta = beta;

    if ((alg == ProcessNoiseAlgorithm::ADMC) || (alg == ProcessNoiseAlgorithm::ASNC)) {
      cfg.use_adaptive_process_noise = true;
    } else {  // Fixed Process Noise parameters
      cfg.use_adaptive_process_noise = false;
    }

    cfg.tau_ambiguity = filter_node["tau_ambiguity"].as<double>(1e6);
    cfg.q_ambiguity = filter_node["q_ambiguity"].as<double>(0.1);

    // Outliers ---------------------------------------------------------------------------------
    cfg.outlier_rejection = filter_node["outlier_rejection"].as<bool>(true);
    cfg.outlier_threshold = filter_node["outlier_threshold"].as<double>(50.0);

    // Covariances ---------------------------------------------------------------------------------
    cfg.sigma_r0_sat = filter_node["sigma_r0_sat"].as<double>(100.0);
    cfg.sigma_v0_sat = filter_node["sigma_v0_sat"].as<double>(1.0);
    cfg.sigma_clk0_sat = filter_node["sigma_clk0_sat"].as<double>(100.0);
    cfg.sigma_clkd0_sat = filter_node["sigma_clkd0_sat"].as<double>(1.0);

    // State Numbers
    // ---------------------------------------------------------------------------------
    cfg.N_t_filter = filter_node["N_t_filter"].as<int>(static_cast<int>(-1));
    if (cfg.N_t_filter == -1) {
      cfg.N_t_filter = static_cast<int>(tspan.size());
    }

    cfg.N_sat_filter = 1;  // single satellite for now
    cfg.N_clk = filter_node["N_clk"].as<int>(2);

    // ------------------------------------------------------------------------------------
    // Add Measurement Flags
    // ---------------------------------------------------------------------------------
    auto meas_config = sim_config.config["measurements"];
    cfg.use_pr = AddMeasFlags(meas_config, "code", MeasType::CODE, cfg.meas_flags);
    cfg.use_carrier = AddMeasFlags(meas_config, "carrier", MeasType::CARRIER, cfg.meas_flags);
    cfg.use_graphic = AddMeasFlags(meas_config, "graphic", MeasType::GRAPHIC, cfg.meas_flags);
    cfg.use_tdcp = AddMeasFlags(meas_config, "tdcp", MeasType::TDCP, cfg.meas_flags);
    cfg.use_ionofree_pr
        = AddMeasFlags(meas_config, "ionofree_code", MeasType::IONOFREE_CODE, cfg.meas_flags);
    cfg.use_ionofree_carrier
        = AddMeasFlags(meas_config, "ionofree_carrier", MeasType::IONOFREE_CARRIER, cfg.meas_flags);
    cfg.use_ionofree_tdcp
        = AddMeasFlags(meas_config, "ionofree_tdcp", MeasType::IONOFREE_TDCP, cfg.meas_flags);
    cfg.use_ionofree_any = cfg.use_ionofree_pr || cfg.use_ionofree_carrier || cfg.use_ionofree_tdcp;
    std::cout << "Use Iono-Free Measurements: " << (cfg.use_ionofree_any ? "True" : "False")
              << std::endl;

    // ------------------------------------------------------------------------------------
    // Ambiguity Measurement Frequencies
    // ------------------------------------------------------------------------------------
    std::vector<MeasFreq> meas_freqs_all = {MeasFreq::L1, MeasFreq::L5, MeasFreq::L1_L5};
    cfg.ambiguity_meas_freqs.clear();

    for (const auto& freq : meas_freqs_all) {
      // Use CARRIER measurements for the frequency -> add to ambiguity measurement frequencies
      if (cfg.meas_flags.find({freq, MeasType::CARRIER}) != cfg.meas_flags.end()) {
        if (cfg.meas_flags[{freq, MeasType::CARRIER}].use_meas) {
          cfg.ambiguity_meas_freqs.push_back(freq);
        }
      }

      // Use GRAPHIC measurements for the frequency -> add to ambiguity measurement frequencies
      if (cfg.meas_flags.find({freq, MeasType::GRAPHIC}) != cfg.meas_flags.end()) {
        if (cfg.meas_flags[{freq, MeasType::GRAPHIC}].use_meas) {
          cfg.ambiguity_meas_freqs.push_back(freq);
        }
      }
      // Use IONOFREE_CARRIER measurements for the frequency -> add to ambiguity measurement
      // frequencies
      if (cfg.meas_flags.find({freq, MeasType::IONOFREE_CARRIER}) != cfg.meas_flags.end()) {
        if (cfg.meas_flags[{freq, MeasType::IONOFREE_CARRIER}].use_meas) {
          cfg.ambiguity_meas_freqs.push_back(freq);
        }
      }
    }
    cfg.n_ambiguity_signal_types = static_cast<int>(cfg.ambiguity_meas_freqs.size());

    // Identifty if we need to estimate float ambiguties ---------------------------------------
    cfg.est_integers = false;
    if ((cfg.use_graphic || cfg.use_carrier || cfg.use_ionofree_carrier)) {
      cfg.est_integers = true;  // We need to estimate float ambiguities for GRAPHIC/carrier
    }

    // ---------------------------------------------------------------------------------
    // Print Ambiguity Measurement Frequencies
    // ---------------------------------------------------------------------------------
    if (verbose) {
      std::cout << "\n[Measurement Flags]\n";
      std::cout << "Number of Ambiguity Signal Types: " << cfg.n_ambiguity_signal_types
                << std::endl;

      // Table header
      std::cout << std::left << std::setw(10) << "Freq" << std::setw(20) << "Meas Type"
                << std::setw(8) << "Use" << std::setw(15) << "Min Alt [km]" << std::setw(15)
                << "Max Alt [km]"
                << "\n";

      std::cout << std::string(68, '-') << "\n";

      for (const auto& [freq_type_pair, flag] : cfg.meas_flags) {
        // Frequency string
        std::string freq_str = (freq_type_pair.first == MeasFreq::L1)   ? "L1"
                               : (freq_type_pair.first == MeasFreq::L5) ? "L5"
                                                                        : "L1+L5";

        // Measurement type string
        std::string meas_type_str
            = (freq_type_pair.second == MeasType::CODE)               ? "CODE"
              : (freq_type_pair.second == MeasType::CARRIER)          ? "CARRIER"
              : (freq_type_pair.second == MeasType::GRAPHIC)          ? "GRAPHIC"
              : (freq_type_pair.second == MeasType::TDCP)             ? "TDCP"
              : (freq_type_pair.second == MeasType::IONOFREE_CODE)    ? "IONOFREE_CODE"
              : (freq_type_pair.second == MeasType::IONOFREE_CARRIER) ? "IONOFREE_CARRIER"
              : (freq_type_pair.second == MeasType::IONOFREE_TDCP)    ? "IONOFREE_TDCP"
                                                                      : "UNDEFINED";

        std::cout << std::left << std::setw(10) << freq_str << std::setw(20) << meas_type_str
                  << std::setw(8) << (flag.use_meas ? "Yes" : "No") << std::setw(15) << std::fixed
                  << std::setprecision(2) << flag.min_alt / 1000 << std::setw(15) << std::fixed
                  << std::setprecision(2) << flag.max_alt / 1000 << "\n";
      }

      std::cout << std::endl;
    }

    // SRP dynamics
    // ---------------------------------------------------------------------------------
    cfg.use_srp = filter_node["use_srp"].as<bool>(false);
    auto dynamics_node = sim_config.config["dynamics"];
    double cr = dynamics_node["cr"].as<double>(1.8);
    double area = dynamics_node["area"].as<double>(1.0);
    double mass = dynamics_node["mass"].as<double>(850.0);
    cfg.true_srp_b_coeff = cr * area / mass;
    cfg.use_srp = filter_node["use_srp"].as<bool>(false);
    cfg.srp_initial_error_ratio = filter_node["srp_initial_error_ratio"].as<double>(0.2);

    // Clock models
    // ---------------------------------------------------------------------------------
    std::string clk_model_str = sim_config.config["clock"].as<std::string>("OCXO");

    // Clock for estimation
    cfg.clk_model_sat = StringToClockModel(clk_model_str);
    cfg.clock_ = std::make_shared<Clock>();
    cfg.clock_->SetModel(cfg.clk_model_sat);
    cfg.clock_dynamics_ = cfg.clock_->GetDynamicsSharedPtr();
    cfg.clock_dynamics_->SetAddNoise(false);  // No noise for estimation clock

    // Clock for true dynamics
    cfg.clock_truth_ = std::make_shared<Clock>();
    cfg.clock_truth_->SetModel(cfg.clk_model_sat);
    cfg.clock_dynamics_truth_ = cfg.clock_truth_->GetDynamicsSharedPtr();
    cfg.clock_dynamics_truth_->SetModel(cfg.clk_model_sat);
    cfg.clock_dynamics_truth_->SetSeed(seed);
    cfg.clock_dynamics_truth_->SetAddNoise(true);  // No noise for true clock

    cfg.inflate_q2 = filter_node["inflate_q2"].as<double>(1.0);

    // Cycle slip ---------------------------------------------------------------------------------
    cfg.cycle_slip_ratio = meas_config["cycle_slip"]["ratio"].as<double>(0.0);
    cfg.cycle_slip_magnitude = meas_config["cycle_slip"]["magnitude"].as<int>(5);
    cfg.cycle_slip_detect_sigma = meas_config["cycle_slip"]["detect_sigma"].as<double>(3.0);
    cfg.cycle_slip_max_cn0 = meas_config["cycle_slip"]["max_cn0"].as<double>(25.0);
    cfg.fix_cycleslip = meas_config["cycle_slip"]["fix_cycleslip"].as<bool>(false);

    // Noise Models
    // ---------------------------------------------------------------------------------
    cfg.sigma_ure = filter_node["sigma_ure"].as<double>(1.0);                        // [m]
    cfg.sigma_iono = filter_node["sigma_iono"].as<double>(1.0);                      // [m]
    cfg.sigma_tdcp = filter_node["sigma_tdcp"].as<double>(0.01);                     // [m]
    cfg.sigma_graphic = filter_node["sigma_graphic"].as<double>(0.1);                // [m]
    cfg.sigma_tdcp_ionofree = filter_node["sigma_tdcp_ionofree"].as<double>(0.001);  // [m]

    cfg.sigma_inflation_code = filter_node["sigma_inflation"]["code"].as<double>(1.0);
    cfg.sigma_inflation_carrier = filter_node["sigma_inflation"]["carrier"].as<double>(1.0);
    cfg.sigma_inflation_tdcp = filter_node["sigma_inflation"]["tdcp"].as<double>(1.0);
    cfg.sigma_inflation_graphic = filter_node["sigma_inflation"]["graphic"].as<double>(1.0);
    cfg.sigma_inflation_ionofree_code
        = filter_node["sigma_inflation"]["ionofree_code"].as<double>(1.0);
    cfg.sigma_inflation_ionofree_carrier
        = filter_node["sigma_inflation"]["ionofree_carrier"].as<double>(1.0);

    // Autodiff ---------------------------------------------------------------------------------
    cfg.use_meas_autodiff = filter_node["use_meas_autodiff"].as<bool>(false);

    // UDU ---------------------------------------------------------------------------------
    cfg.use_udu = filter_node["use_udu"].as<bool>(false);
    cfg.udu_predict_single = filter_node["udu_predict_single"].as<bool>(false);

    // Smoothing ---------------------------------------------------------------------------------
    cfg.N_smooth_iter = sim_config.N_smooth_iter;
    cfg.smooth_fraction = filter_node["smoothing_fraction"].as<double>(0.33);
    cfg.run_smoother = filter_node["run_smoother"].as<bool>(false);
    cfg.smooth_start_tidx = int(cfg.N_t_filter * (1 - cfg.smooth_fraction));

    // Dynamics
    // ---------------------------------------------------------------------------------
    cfg.joint_orbit_clock_dynamics = filter_node["joint_orbit_clock_dynamics"].as<bool>(true);

    // Recompute
    cfg.recompute = sim_config.recompute_filter;

    // verbose settings validation
    // ---------------------------------------------------------------------------------
    auto verbose_node = sim_config.config["verbose"];
    cfg.verbose_filter_dynamics = verbose_node["filter_dynamics"].as<int>(0);
    cfg.verbose_filter_predict = verbose_node["filter_predict"].as<int>(0);
    cfg.verbose_filter_measurement = verbose_node["filter_measurement"].as<int>(1000);
    cfg.verbose_filter_update = verbose_node["filter_update"].as<int>(1000);
    cfg.verbose_filter_smoothing = verbose_node["filter_smoothing"].as<int>(1000);

    if (cfg.cycle_slip_ratio < 0.0 || cfg.cycle_slip_ratio > 1.0) {
      throw std::runtime_error("Cycle slip ratio must be between 0 and 1");
    }

    std::cout << "Filter Configuration Finished" << std::endl;

    return cfg;
  }

  std::string ClockModelToString(ClockModel model) {
    switch (model) {
      case ClockModel::OCXO: return "OCXO";
      case ClockModel::USO: return "USO";
      case ClockModel::CSAC: return "CSAC";
      case ClockModel::MINI_RAFS: return "MINI_RAFS";
      case ClockModel::RAFS: return "RAFS";
      case ClockModel::DSAC: return "DSAC";
      default: return "UNDEFINED";
    }
  }

  FilterState InitializeFilterState(const FilterConfig& cfg, const SimulationConfig& sim_config,
                                    MatXd sat_rv_mci, VecXd tspan, double t0_tai, int seed,
                                    bool verbose) {
    FilterState state;

    state.N_rva = cfg.dmc ? 9 : 6;
    state.N_srp = cfg.use_srp ? 1 : 0;
    state.N_clk = cfg.N_clk;
    state.N_ambiguity_est = 0;   // temporary, will be set later based on measurements
    state.N_ambiguity_full = 1;  // temporary, will be set later based on measurements
    state.N_ekf_sat_est = state.N_rva + state.N_clk + state.N_srp + state.N_ambiguity_est;
    state.N_t = cfg.N_t_filter;

    // Initialize state arrays
    state.t0_tai = t0_tai;
    state.ts_filter = tspan.head(cfg.N_t_filter);
    state.rva_true_sat.resize(cfg.N_t_filter, state.N_rva);
    state.rva_est_sat.resize(cfg.N_t_filter, state.N_rva);
    state.rva_sigma_sat.resize(cfg.N_t_filter, state.N_rva);
    state.clk_true_sat.resize(cfg.N_t_filter, cfg.N_clk);
    state.clk_est_sat.resize(cfg.N_t_filter, cfg.N_clk);
    state.clk_sigma_sat.resize(cfg.N_t_filter, cfg.N_clk);
    state.srp_true_sat.resize(cfg.N_t_filter);
    state.srp_est_sat.resize(cfg.N_t_filter);
    state.srp_sigma_sat.resize(cfg.N_t_filter);
    state.ambiguity_true_sat.resize(cfg.N_t_filter, state.N_ambiguity_full);
    state.ambiguity_est_sat.resize(cfg.N_t_filter, state.N_ambiguity_full);
    state.ambiguity_sigma_sat.resize(cfg.N_t_filter, state.N_ambiguity_full);
    state.ambiguity_is_estimated.resize(cfg.N_t_filter, state.N_ambiguity_full);
    state.P_est = MatXd::Zero(state.N_ekf_sat_est, state.N_ekf_sat_est);
    // RTN errors
    state.rtn_err_pos.resize(cfg.N_t_filter, 3);
    state.rtn_err_vel.resize(cfg.N_t_filter, 3);
    state.rtn_err_pos_sigma.resize(cfg.N_t_filter, 3);
    state.rtn_err_vel_sigma.resize(cfg.N_t_filter, 3);
    // Measurements
    state.num_meas.resize(cfg.N_t_filter);

    // Setup clock dynamics -------------------------------------------------
    RandomEngine::SetSeed(seed);

    // Initialize satellite states
    VecXd sigmas_rva0 = VecX::Zero(state.N_rva);
    if (state.N_rva == 6) {
      sigmas_rva0 << cfg.sigma_r0_sat * Vec3::Ones(), cfg.sigma_v0_sat * Vec3::Ones();
    } else {
      sigmas_rva0 << cfg.sigma_r0_sat * Vec3::Ones(), cfg.sigma_v0_sat * Vec3::Ones(),
          cfg.sigma_v0_sat / 1e3 * Vec3::Ones();
    }
    MatXd cov_rva = sigmas_rva0.array().square().matrix().asDiagonal();

    // First 6 columns: position and velocity
    state.rva_true_sat.leftCols(6) = sat_rv_mci.topRows(cfg.N_t_filter);

    // state.rva_est_sat.row(0).head(6)
    state.rva_est_sat.row(0).head(6)
        = SampleMvNormal(state.rva_true_sat.row(0).head(6), cov_rva, 1).transpose();
    state.rva_sigma_sat.row(0) = sigmas_rva0;
    state.P_est.block(0, 0, state.N_rva, state.N_rva) = cov_rva;

    // Initialize satellite clock ---------------------------------------------------------
    VecXd sigmas_clk0 = VecX::Zero(cfg.N_clk);
    sigmas_clk0(0) = cfg.sigma_clk0_sat;
    sigmas_clk0(1) = cfg.sigma_clkd0_sat;
    if (cfg.N_clk == 3) sigmas_clk0(2) = cfg.sigma_clkd0_sat / 1e3;
    MatXd cov_clk = sigmas_clk0.array().square().matrix().asDiagonal();

    // Propagate true clock states ------------------------------------------------
    if (cfg.joint_orbit_clock_dynamics) {
      // In this version, we include relativistic clock drift in the clock propagation
      state.clk_true_sat.row(0) = VecX::Zero(cfg.N_clk);  // start at zero clock offset/rate

      // Load data from cache
      std::string orbit_config_str = "norbit_"
                                     + std::to_string(static_cast<int>(sim_config.N_orbit)) + "_dt"
                                     + std::to_string(static_cast<int>(sim_config.dt)) + "s_sat"
                                     + std::to_string(sim_config.sat_id);
      std::string clock_filename = "true_clock_" + ClockModelToString(cfg.clk_model_sat) + "_seed_"
                                   + std::to_string(seed) + ".h5";
      std::filesystem::path output_path
          = GetOutputDir("iono_delay") / "clocks" / orbit_config_str / clock_filename;
      H5Easy::File clock_file(output_path.string(), H5Easy::File::ReadOnly);

      auto func = []() { return MatX3(); };
      std::cout << "Loading true clock data from: " << output_path << std::endl;
      state.clk_true_sat = LoadOrRecompute<-1, 3, double>("/clk_true_sat", clock_file, false, func);
      std::cout << "True clock data loaded. Size: " << state.clk_true_sat.rows() << " x "
                << state.clk_true_sat.cols() << std::endl;

    } else {
      state.clk_true_sat
          = cfg.clock_dynamics_truth_->Propagate(VecX::Zero(cfg.N_clk), state.ts_filter).array()
            * C;
    }

    state.clk_est_sat.row(0) = SampleMvNormal(state.clk_true_sat.row(0), cov_clk, 1).transpose();
    state.clk_sigma_sat.row(0) = sigmas_clk0;
    state.P_est.block(state.N_rva, state.N_rva, cfg.N_clk, cfg.N_clk) = cov_clk;

    // Initialize SRP coefficient if used
    if (cfg.use_srp) {
      state.srp_true_sat.setConstant(cfg.true_srp_b_coeff);
      double sigma_srp = cfg.true_srp_b_coeff * cfg.srp_initial_error_ratio;
      state.srp_est_sat(0) = SampleNormal(cfg.true_srp_b_coeff, sigma_srp);
      state.srp_sigma_sat(0) = sigma_srp;
      state.P_est(state.N_rva + cfg.N_clk, state.N_rva + cfg.N_clk) = sigma_srp * sigma_srp;
    }

    // Initialize satellite ambiguities to zero for now
    state.ambiguity_true_sat.setZero();
    state.ambiguity_est_sat.setZero();
    state.ambiguity_sigma_sat.setZero();
    state.ambiguity_is_estimated.setZero();

    if (verbose) {
      std::cout << " " << std::endl;
      std::cout << "[Filter State Initialization Summary]" << std::endl;
      std::cout << "-- RVa State --" << std::endl;
      std::cout << "  True RVa:\n" << state.rva_true_sat.row(0).head(6) << std::endl;
      std::cout << "  Estimated RVa:\n" << state.rva_est_sat.row(0).head(6) << std::endl;
      std::cout << "  Est Err RVa:\n"
                << (state.rva_est_sat.row(0).head(6) - state.rva_true_sat.row(0).head(6))
                << std::endl;
      std::cout << "-- Clock State --" << std::endl;
      std::cout << "  True Clock:\n" << state.clk_true_sat.row(0) << std::endl;
      std::cout << "  Estimated Clock:\n" << state.clk_est_sat.row(0) << std::endl;
      std::cout << "  Est Err Clock:\n"
                << (state.clk_est_sat.row(0) - state.clk_true_sat.row(0)) << std::endl;
      if (cfg.use_srp) {
        std::cout << "-- SRP Coefficient --" << std::endl;
        std::cout << "  True SRP Coeff: " << state.srp_true_sat(0) << " m^2/kg" << std::endl;
        std::cout << "  Estimated SRP Coeff: " << state.srp_est_sat(0) << " m^2/kg" << std::endl;
        std::cout << "  Est Err SRP Coeff: " << (state.srp_est_sat(0) - state.srp_true_sat(0))
                  << " m^2/kg" << std::endl;
        std::cout << " " << std::endl;
      }
    }

    return state;
  }

  Ptr<NBodyDynamics> SetupDynamics(bool use_ad, bool use_srp, int sph_order, double srp_coeff,
                                   double dt_prop) {
    // Setup dynamics
    Ptr<NBodyDynamics> dynamics = std::make_shared<NBodyDynamics>();
    dynamics->SetIntegrator(IntegratorType::RKF45);
    dynamics->SetIntegrator(IntegratorType::RKF45);
    dynamics->AddBody(Body::Moon(sph_order, sph_order));
    dynamics->AddBody(Body::Earth());
    dynamics->AddBody(Body::Sun());
    dynamics->SetFrame(Frame::MOON_CI);
    dynamics->SetTimeStep(dt_prop);
    dynamics->SetAutodiff(true);
    if (use_srp) {
      dynamics->SetSrpCoeff(srp_coeff);
    }

    return dynamics;
  }

  JointState CreateJointState(const SimulationConfig& sim_config, const FilterConfig& filter_cfg,
                              FilterState& state, bool verbose) {
    // Setup dynamics
    SetLupntEpoch(state.t0_tai);
    Ptr<NBodyDynamics> dyn_filter
        = SetupDynamics(true, filter_cfg.use_srp, 18, state.srp_est_sat(0), sim_config.dt_prop);

    // Test Propagation
    bool test_propagation = false;

    if (test_propagation) {
      std::cout << " " << std::endl;
      std::cout << "  [Propagation Test]" << std::endl;
      int N_prop = 10000;
      Cart6 rv0_true_ci = Cart6(state.rva_true_sat.row(0).transpose(), Frame::MOON_CI);
      Cart6 rvf_true_ci = Cart6(state.rva_true_sat.row(N_prop).transpose(), Frame::MOON_CI);

      Cart6 rvf_propagated_ci
          = dyn_filter->Propagate(rv0_true_ci, 0, state.ts_filter[N_prop], nullptr);

      std::cout << "  Propagation Time: " << state.ts_filter[N_prop] - state.ts_filter[0] << " s"
                << std::endl;
      std::cout << "  True Final State (CI):\n" << rvf_true_ci << std::endl;
      std::cout << "  Propagated Final State (CI):\n" << rvf_propagated_ci << std::endl;
      std::cout << "  Estimation Error (CI):\n" << rvf_propagated_ci - rvf_true_ci << std::endl;
      std::cout << " " << std::endl;
    }

    // Setup joint state
    Cart6 rv0_ci = Cart6(state.rva_est_sat.row(0).transpose(), Frame::MOON_CI);
    // Parameter State
    ParamState param_state = dyn_filter->GetParams();

    // Clock state
    if (filter_cfg.N_clk == 3) {
      ClockState3 clock_state3;
      clock_state3.b() = state.clk_est_sat(0, 0);   // [s]
      clock_state3.d() = state.clk_est_sat(0, 1);   // [s/s]
      clock_state3.dr() = state.clk_est_sat(0, 2);  // [s/s^2]
      filter_cfg.clock_->SetState(clock_state3);
    } else {
      ClockState2 clock_state2;
      clock_state2.b() = state.clk_est_sat(0, 0);  // [s]
      clock_state2.d() = state.clk_est_sat(0, 1);  // [s/s]
      filter_cfg.clock_->SetState(clock_state2);
    }

    // Process Noise
    MatXd Q_a = filter_cfg.Q_a;        // [m/s^2/sqrt(Hz)]
    VecXd beta = filter_cfg.dmc_beta;  // For ADMC

    auto process_noise_rv = [Q_a, beta](const State& x, Real t0, Real tf) -> MatXd {
      int n = x.size();
      MatXd Q = MatXd::Zero(n, n);
      Real dt = tf - t0;
      MatXd Q_rv;

      if (n == 6) {
        Q_rv = ProcessNoisePosVel(Q_a, dt);
      } else if (n == 9) {
        Q_rv = ProcessNoisePosVelAcc(Q_a, dt, beta);
      }

      return Q_rv;
    };

    ClockModel clk_model = filter_cfg.clk_model_sat;
    Ptr<ClockDynamics> clock_dynamics = filter_cfg.clock_dynamics_;
    int N_clk = filter_cfg.N_clk;
    auto proc_noise_clock
        = [clk_model, clock_dynamics, N_clk](const State& x, Real t0, Real tf) -> MatXd {
      if (N_clk == 3) {
        return clock_dynamics->ThreeStateNoise(clk_model, tf - t0);
      } else {
        return clock_dynamics->TwoStateNoise(clk_model, tf - t0);
      }
    };

    auto proc_noise_rv_ptr = std::make_shared<ProcessNoiseFunction>(process_noise_rv);
    auto proc_noise_clock_ptr = std::make_shared<ProcessNoiseFunction>(proc_noise_clock);

    // Define Joint State ------------------------------------------------
    JointState joint_state;
    if (filter_cfg.joint_orbit_clock_dynamics) {
      int N_rva = state.N_rva;
      auto proc_noise_rvc = [process_noise_rv, proc_noise_clock, N_rva, N_clk](
                                const State& x, Real t0, Real tf) -> MatXd {
        MatXd Q_rv = process_noise_rv(x.head(N_rva), t0, tf);
        MatXd Q_clk = proc_noise_clock(x.tail(N_clk), t0, tf);
        MatXd Q = MatXd::Zero(N_rva + N_clk, N_rva + N_clk);
        Q.topLeftCorner(N_rva, N_rva) = Q_rv;
        Q.bottomRightCorner(N_clk, N_clk) = Q_clk;
        return Q;
      };
      auto proc_noise_rvc_ptr = std::make_shared<ProcessNoiseFunction>(proc_noise_rvc);

      JointOrbitClockState x_joint = JointOrbitClockState(rv0_ci, filter_cfg.clock_->GetState());
      Ptr<JointOrbitClockDynamics> joint_dynamics = std::make_shared<JointOrbitClockDynamics>();
      joint_dynamics->SetOrbitDynamics(dyn_filter);
      joint_dynamics->SetClockDynamics(clock_dynamics);
      joint_dynamics->SetAddNoise(false);
      joint_state.Add(x_joint, joint_dynamics, std::move(proc_noise_rvc_ptr), param_state,
                      std::vector<EstType>{ESTIMATED, FIXED});
    } else {
      if (filter_cfg.use_srp) {
        // When using SRP, estimate them
        joint_state.Add(rv0_ci, dyn_filter, std::move(proc_noise_rv_ptr), param_state,
                        std::vector<EstType>{ESTIMATED, FIXED});
      } else {
        joint_state.Add(rv0_ci, dyn_filter, std::move(proc_noise_rv_ptr), param_state,
                        std::vector<EstType>{FIXED, FIXED});
      }
      // Do not add clock state to joint state for now, handle clock separately
      joint_state.Add(filter_cfg.clock_->GetState(), clock_dynamics,
                      std::move(proc_noise_clock_ptr), ParamState(0), std::vector<EstType>{});
    }

    // Summary ------------------------------------------------------
    if (verbose) {
      std::cout << "  " << std::endl;
      std::cout << "[Joint State Summary]" << std::endl;
      std::cout << "  Joint State STM   Size: " << joint_state.GetStmSize() << std::endl;
      std::cout << "  Joint State State Size: " << joint_state.GetStateSize() << std::endl;
      std::cout << "  Joint State Param Size: " << joint_state.GetParamSize() << std::endl;
      std::cout << "     - Estimated Params  : " << joint_state.GetEstimatedParamSize()
                << std::endl;
      std::cout << "     - Considered Params : " << joint_state.GetConsideredParamSize()
                << std::endl;
      std::cout << "     - Fixed Params      : " << joint_state.GetFixedParamSize() << std::endl;
      std::cout << "-- State Vector --" << std::endl;
      std::cout << joint_state.GetState().transpose() << std::endl;
      std::cout << "-- Parameter Vector --" << std::endl;
      std::cout << joint_state.GetParams().transpose() << std::endl;
      std::cout << "-- STM State Vector --" << std::endl;
      std::cout << joint_state.GetStmState().transpose() << std::endl;
      std::cout << " " << std::endl;
    }

    return joint_state;
  }

}  // namespace filtering_sim
