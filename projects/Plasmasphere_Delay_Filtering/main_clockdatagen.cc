
#include <lupnt/lupnt.h>

#include "src/dynamics_models.h"
#include "src/filter_setup.h"

using namespace lupnt;
using std::cout;
using std::endl;
using namespace filtering_sim;

struct OrbitData {
  VecXd tspan;
  VecXd tais;
  MatX6d posvel_rx_mci;
  MatX6d posvel_rx_gcrf;
  MatX6d posvel_rx_pa;
  int n_steps;
  std::string orbit_config_str;
};

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

OrbitData load_orbit(std::filesystem::path cache_path) {
  bool recompute = false;

  H5Easy::File rx_cache_file(cache_path.string(), H5Easy::File::ReadOnly);

  std::cout << "Loading receiver history from cache..." << std::endl;

  // Dummy lambda for compatibility
  auto func = []() { return MatX6(); };
  auto func_vec = []() { return VecX(); };

  OrbitData orbit_data;

  orbit_data.tspan = LoadOrRecompute<-1, 1, double>("/tspan", rx_cache_file, recompute, func_vec);
  orbit_data.tais = LoadOrRecompute<-1, 1, double>("/t_tai", rx_cache_file, recompute, func_vec);
  orbit_data.posvel_rx_mci
      = LoadOrRecompute<-1, 6, double>("/posvel_rx_mci", rx_cache_file, recompute, func);
  orbit_data.posvel_rx_gcrf
      = LoadOrRecompute<-1, 6, double>("/posvel_rx_gcrf", rx_cache_file, recompute, func);
  orbit_data.posvel_rx_pa
      = LoadOrRecompute<-1, 6, double>("/posvel_rx_pa", rx_cache_file, recompute, func);

  orbit_data.n_steps = static_cast<int>(orbit_data.tspan.size());

  std::cout << "  Loaded " << orbit_data.n_steps << " orbit steps from cache." << std::endl;

  return orbit_data;
};

void PropagateTrueClockState(OrbitData orbit_data, ClockModel clk_model, int seed,
                             bool verbose = true) {
  // Output Path
  std::string orbit_config_str = orbit_data.orbit_config_str;
  std::string output_filename
      = "true_clock_" + ClockModelToString(clk_model) + "_seed_" + std::to_string(seed) + ".h5";
  std::filesystem::path output_path
      = GetOutputDir("iono_delay") / "clocks" / orbit_config_str / output_filename;

  // Skip if already exists
  if (std::filesystem::exists(output_path)) {
    std::cout << " True clock state file already exists at: " << output_path << std::endl;
    std::cout << " Skipping propagation..." << std::endl;
    return;
  }

  // Orbit Data
  VecXd ts_filter = orbit_data.tspan;
  int n_steps = orbit_data.n_steps;
  MatX6d rva_true_sat = orbit_data.posvel_rx_mci;

  int N_t = rva_true_sat.rows();
  Cart6 rv0_ci_true = Cart6(rva_true_sat.row(0).transpose(), Frame::MOON_CI);

  double cr = 1.8;
  double area = 1.0;
  double mass = 850.0;
  double srp_coeff = cr * area / mass;

  int N_rva = 6;        // RVa State Size
  int N_clk = 3;        // Filter Clock State
  bool use_srp = true;  // Whether to use SRP in dynamics

  // Clock Dynamics -------------------------------------------------------------
  Ptr<Clock> clock_truth = std::make_shared<Clock>();
  clock_truth->SetModel(clk_model);
  Ptr<ClockDynamics> clock_dynamics_truth = clock_truth->GetDynamicsSharedPtr();
  clock_dynamics_truth->SetModel(clk_model);
  clock_dynamics_truth->SetSeed(seed);
  clock_dynamics_truth->SetAddNoise(true);  // noise for true clock

  State true_clock_state;
  if (N_clk == 3) {
    ClockState3 clk_state3;
    clk_state3.b() = 0.0;   // [s]
    clk_state3.d() = 0.0;   // [s/s]
    clk_state3.dr() = 0.0;  // [s/s^2]
    true_clock_state = clk_state3;
  } else {
    ClockState2 clk_state2;
    clk_state2.b() = 0.0;  // [s]
    clk_state2.d() = 0.0;  // [s/s]
    true_clock_state = clk_state2;
  }

  // Dynamics
  double dt_prop = 1.0;  // Propagation time step [s]
  Ptr<NBodyDynamics> dyn_truth = SetupDynamics(false, use_srp, 50, srp_coeff, dt_prop);

  // Propagate the Clock and State to recompute true clock state with relativistic effects
  JointOrbitClockState x_joint_true = JointOrbitClockState(rv0_ci_true, true_clock_state);
  Ptr<JointOrbitClockDynamics> joint_dynamics_truth = std::make_shared<JointOrbitClockDynamics>();
  joint_dynamics_truth->SetOrbitDynamics(dyn_truth);
  joint_dynamics_truth->SetClockDynamics(clock_dynamics_truth);
  joint_dynamics_truth->SetAddNoise(true);

  std::cout << " Propagating True Clock State with Relativistic Effects..." << std::endl;
  MatX x_joint_propagated = joint_dynamics_truth->Propagate(x_joint_true, ts_filter, nullptr);
  MatXd clk_true_sat = x_joint_propagated.block(0, N_rva, N_t, N_clk);

  // Check if final propagated orbit matches the stored true orbit
  if (verbose) {
    VecX rvf_ci = x_joint_propagated.row(x_joint_propagated.rows() - 1).head(6);
    VecX rvf_ci_true_check = rva_true_sat.row(N_t - 1).transpose();
    VecX err_rvf_ci_check = rvf_ci - rvf_ci_true_check;

    auto old_format = std::cout.flags();
    auto old_precision = std::cout.precision();
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  [True Clock State Propagation]" << std::endl;
    std::cout << "    True Final State Error (CI):\n" << err_rvf_ci_check.transpose() << std::endl;
    std::cout << "    Final Clock Bias (Propagated) :\n"
              << x_joint_propagated.row(N_t - 1).segment(6, N_clk) << std::endl;
    std::cout << "    Propagated Clock State Size: " << clk_true_sat.rows() << " x "
              << clk_true_sat.cols() << std::endl;
    std::cout << " " << std::endl;
    std::cout.flags(old_format);
    std::cout.precision(old_precision);
  }

  // Update the clock state in the filter config

  // Save to h5 file
  // Generate Directory if not exists
  std::filesystem::create_directories(output_path.parent_path());

  // Save Results
  H5Easy::File results_file = GetH5File(output_path, true);
  Dump(results_file, "/clk_true_sat", clk_true_sat);

  std::cout << " True clock state saved to: " << output_path << std::endl;
}

int main(int argc, char* argv[]) {
  // Settings
  int num_mc = 50;
  int sat_id = 0;
  int N_orbit = 6;
  int dt = 1;  // seconds
  ClockModel clock_model = ClockModel::OCXO;

  // Orbit Data
  std::string rx_filename = "lcrns_sat" + std::to_string(sat_id) + "_norbit"
                            + std::to_string(static_cast<int>(N_orbit)) + "_dt"
                            + std::to_string(static_cast<int>(dt)) + ".h5";
  std::filesystem::path rx_cache_path = GetOutputDir("iono_delay") / "orbits" / "h5" / rx_filename;
  OrbitData orbit_data = load_orbit(rx_cache_path);

  orbit_data.orbit_config_str = "norbit_" + std::to_string(static_cast<int>(N_orbit)) + "_dt"
                                + std::to_string(static_cast<int>(dt)) + "s_sat"
                                + std::to_string(sat_id);

  SetLupntEpoch(orbit_data.tais(0));

  // Run Monte Carlo Simulations
  auto pbar = Logger::GetProgressBar(num_mc, "Running MC", "Clock MC");

  for (int i = 0; i < num_mc; i++) {
    std::cout << "======================" << std::endl;
    std::cout << " Monte Carlo Run " << i + 1 << " / " << num_mc << std::endl;
    std::cout << "======================" << std::endl;

    // Set random seed for each MC run
    int seed = i;

    // Execute filtering simulation
    PropagateTrueClockState(orbit_data, clock_model, seed, true);

    pbar->Update();
  }
  pbar->Finish();
}
