/**
 * @file example_ekf_elfo_gps.cc
 * @author Stanford NAV Lab
 * @brief   Example of EKF for lunar orbit estimation
 * @version 0.1
 * @date 2024-01-19
 *
 * @copyright Copyright (c) 2024
 *
 */

// lupnt includes
#include <lupnt/lupnt.h>

#include <vector>

using namespace lupnt;
namespace sp = spice;

class MyFile {
private:
  std::ofstream file;

public:
  MyFile(const std::string& filepath) {
    file.open(filepath);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file: " + filepath);
    }
  }

  ~MyFile() {
    if (file.is_open()) {
      file.close();
    }
  }

  // Operator <<
  template <typename T> MyFile& operator<<(const T& data) {
    file << data;
    return *this;
  }
};

template <typename T> class Timestamped {
public:
  Timestamped(double timestamp, const T& data) : timestamp(timestamp), data(data) {}

  double GetTimestamp() const { return timestamp; }
  const T& GetData() const { return data; }

private:
  double timestamp;
  T data;
};

class DataHistory {
public:
  template <typename VectorType>
  void AddData(const std::string& key, double timestamp, const VectorType& data) {
    if (historyData.find(key) == historyData.end()) {
      historyData[key] = {};
    }
    historyData[key].push_back(Timestamped<VecXd>(timestamp, data.template cast<double>()));
  }

  const std::vector<Timestamped<VecXd>>& GetData(const std::string& key) const {
    auto it = historyData.find(key);
    if (it != historyData.end()) {
      return it->second;
    }
    throw std::runtime_error("No data found for the given key.");
  }

  void AddHeader(const std::string& key, const std::string& header) {
    headers[key].push_back(header);
  }

  const std::map<std::string, std::vector<Timestamped<VecXd>>>& GetData() const {
    return historyData;
  }

  const std::map<std::string, std::vector<std::string>>& GetHeaders() const { return headers; }

private:
  std::map<std::string, std::vector<Timestamped<VecXd>>> historyData;
  std::map<std::string, std::vector<std::string>> headers;
};

class FileWriter {
public:
  FileWriter(const std::filesystem::path& basePath, const bool make_dirs = false)
      : basePath(basePath) {
    // Create the base path if it does not exist
    if (make_dirs) {
      std::filesystem::create_directories(basePath);
    }
    // Check if the base path exists
    if (!std::filesystem::exists(basePath)) {
      throw std::runtime_error("Base path does not exist: " + basePath.string());
    }
    // Remove all csv files in the directory
    for (const auto& entry : std::filesystem::directory_iterator(basePath)) {
      if (entry.path().extension() == ".csv") {
        std::filesystem::remove(entry.path());
      }
    }
  }

  void WriteData(const DataHistory& dataHistory) {
    for (const auto& [key, data] : dataHistory.GetData()) {
      std::filesystem::path filepath = basePath / (key + ".csv");

      MyFile file(filepath.string());
      // Check headers
      auto it = dataHistory.GetHeaders().find(key);
      if (it != dataHistory.GetHeaders().end()) {
        file << "t";
        for (const auto& header : it->second) {
          file << "," << header;
        }
        file << "\n";
      }
      for (const auto& timestampedData : data) {
        file << timestampedData.GetTimestamp();
        if (timestampedData.GetData().size() > 0) {
          file << "," << timestampedData.GetData().transpose().format(fmt);
        }
        file << "\n";
      }
    }
  }

private:
  std::filesystem::path basePath;
  Eigen::IOFormat fmt{Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n"};
};

// Util Functions
MatXd ConstructInitCovariance(double pos_err, double vel_err, double clk_bias_err,
                              double clk_drift_err) {
  Mat6d P_rv = Mat6d::Zero();
  P_rv.block(0, 0, 3, 3) = Mat3d::Identity() * pow(pos_err, 2);
  P_rv.block(3, 3, 3, 3) = Mat3d::Identity() * pow(vel_err, 2);

  Mat2d P_clk = Mat2d::Zero();
  P_clk(0, 0) = pow(clk_bias_err, 2);
  P_clk(1, 1) = pow(clk_drift_err, 2);

  MatXd P0 = BlkDiagD(P_rv, P_clk);

  return P0;
};

void AddStateEstimationData(const std::shared_ptr<DataHistory> data_history,
                            const std::shared_ptr<Spacecraft> sat, EKF* ekf,
                            GnssConstellation* gps_const, GnssMeasurement* meas, double t,
                            double epoch) {
  // Navigation
  data_history->AddData("z_true", t, ekf->z_true_);
  data_history->AddData("z_pred", t, ekf->z_pred_);
  data_history->AddData("CN0", t, meas->GetCN0());

  data_history->AddData("vis_earth", t, meas->GetEarthOccultation());
  data_history->AddData("vis_moon", t, meas->GetMoonOccultation());
  data_history->AddData("vis_antenna", t, meas->GetMoonOccultation());
  data_history->AddData("vis_ionos", t, meas->GetMoonOccultation());

  // Moon spacecraft
  CartesianOrbitState state = sat->GetCartesianGCRFStateAtEpoch(epoch);
  CartesianOrbitState state_mi = ConvertOrbitStateFrame(state, epoch, Frame::MOON_CI);
  CartesianOrbitState state_gcrf = ConvertOrbitStateFrame(state, epoch, Frame::GCRF);
  data_history->AddData("rv_moon_mi", t, state_mi.GetVec());
  data_history->AddData("rv_moon_gcrf", t, state_gcrf.GetVec());

  // Estimation
  data_history->AddData("rv", t, sat->GetOrbitState()->GetVec());
  data_history->AddData("rv_pred", t, ekf->xbar_.head(6));
  data_history->AddData("rv_est", t, ekf->x_.head(6));

  data_history->AddData("clk", t, sat->GetClockState().GetVec());
  data_history->AddData("clk_pred", t, ekf->xbar_.tail(2));
  data_history->AddData("clk_est", t, ekf->x_.tail(2));

  data_history->AddData("P_rv", t, ekf->P_.diagonal().segment(0, 6));
  data_history->AddData("P_clk", t, ekf->P_.diagonal().segment(6, 2));

  // GPS constellation
  for (int i = 0; i < gps_const->GetNumSatellites(); i++) {
    CartesianOrbitState sate = gps_const->GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
    CartesianOrbitState sate_mi = ConvertOrbitStateFrame(sate, epoch, Frame::MOON_CI);
    CartesianOrbitState state_gcrf = ConvertOrbitStateFrame(sate, epoch, Frame::GCRF);

    std::string name = "sat" + std::to_string(i);
    data_history->AddData(name + "_mi", t, state.GetVec());
    data_history->AddData(name + "_gcrf", t, state_gcrf.GetVec());
  }

  // Bodies
  Vec6 vz6;
  vz6.setZero();
  data_history->AddData("earth_mi", t, ConvertFrame(epoch, vz6, Frame::GCRF, Frame::MOON_CI));
  data_history->AddData("moon_gcrf", t, ConvertFrame(epoch, vz6, Frame::MOON_CI, Frame::GCRF));
};

void PrintProgressHeader() {
  std::cout << "Run Simulation" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Time [min]  | Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms]" << std::endl;
  std::cout << "--------------------------------------------------------------" << std::endl;
}

VecXd ComputeEstimationErrors(const Ptr<Spacecraft> sat, EKF* ekf) {
  auto x_est = ekf->x_;
  auto x_true = sat->GetStateVec();

  double x_pos_err = 1000 * (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
  double x_vel_err = 1e6 * (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
  double x_clk_bias_err = 3e8 * abs((x_true(6) - x_est(6)).val());
  double x_clk_drift_err = 3e8 * abs((x_true(7) - x_est(7)).val());

  VecXd est_err(4);
  est_err << x_pos_err, x_vel_err, x_clk_bias_err, x_clk_drift_err;

  return est_err;
}

void PrintProgress(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err) {
  std::cout.precision(5);
  std::cout << std::left << std::setw(12) << t / 60 << " " << std::left << std::setw(12)
            << x_pos_err << "  " << std::left << std::setw(14) << x_vel_err << "   " << std::left
            << std::setw(16) << x_clk_bias_err << std::endl;
};

void PrintEKFDebugInfo(int tidx, const Ptr<Spacecraft> sat, EKF* ekf, bool error_only = false) {
  auto x_bar = ekf->xbar_;
  auto x_est = ekf->x_;
  auto x_true = sat->GetStateVec();

  double x_pos_err_bar = 1000 * (x_true.segment(0, 3) - x_bar.segment(0, 3)).norm().val();
  double x_vel_err_bar = 1e6 * (x_true.segment(3, 3) - x_bar.segment(3, 3)).norm().val();
  double x_clk_bias_err_bar = 3e8 * abs((x_true(6) - x_bar(6)).val());

  double x_pos_err = 1000 * (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
  double x_vel_err = 1e6 * (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
  double x_clk_bias_err = 3e8 * abs((x_true(6) - x_est(6)).val());

  std::cout << " ------------------- " << std::endl;
  std::cout << "  Time: " << tidx << std::endl;
  if (!error_only) {
    std::cout << "  Q:  " << std::endl << ekf->Q_ << std::endl;
    std::cout << "  Kalman Gain: " << std::endl << ekf->K_ << std::endl;
    std::cout << "  H:  " << std::endl << ekf->H_ << std::endl;
    std::cout << "  R:  " << std::endl << ekf->R_.diagonal().transpose() << std::endl;
    std::cout << "  S:  " << std::endl << ekf->S_.diagonal().transpose() << std::endl;
    std::cout << " " << std::endl;
  }
  std::cout << "  Meas   Residuals: " << 1000 * ekf->dy_.transpose() << std::endl;
  std::cout << "  Linear Residuals: " << 1000 * (ekf->H_ * ekf->dx_).transpose() << std::endl;
  std::cout << "  dx: " << ekf->dx_.transpose() << std::endl;
  std::cout << "  Pos Err (Bar): " << x_pos_err_bar
            << "  (3sigma : " << 3 * 1000 * sqrt(ekf->Pbar_.block(0, 0, 3, 3).diagonal().trace())
            << " )    Pos Err (Est): " << x_pos_err
            << "  (3sigma: " << 3 * 1000 * sqrt(ekf->P_.block(0, 0, 3, 3).diagonal().trace())
            << " )" << std::endl;
  std::cout << "  Vel Err (Bar): " << x_vel_err_bar
            << "  (3sigma : " << 3 * 1e6 * sqrt(ekf->Pbar_.block(3, 3, 3, 3).diagonal().trace())
            << "  )   Vel Err (Est): " << x_vel_err
            << "  (3sigma: " << 3 * 1e6 * sqrt(ekf->P_.block(3, 3, 3, 3).diagonal().trace()) << " )"
            << std::endl;
  std::cout << "  Clk Err (Bar): " << x_clk_bias_err_bar
            << "  (3sigma : " << 3 * 3e8 * sqrt(ekf->Pbar_.block(6, 6, 1, 1).diagonal().trace())
            << "  )   Clk Err (Est): " << x_clk_bias_err
            << "  (3sigma: " << 3 * 3e8 * sqrt(ekf->P_.block(6, 6, 1, 1).diagonal().trace()) << " )"
            << std::endl;
  std::cout << " ------------------- " << std::endl;
  std::cout << "  " << std::endl;
}

/**
 * @brief Print Estimation Errors
 *
 * @param num_meas (n_time,)   Number of GPS measurements
 * @param error_mat (4, n_time)  Error Mat (Position, Velocity, Clock Bias,
 * Clock Drift)
 */
void PrintEstimationStatistics(VecXd num_meas, MatXd error_mat, double data_ratio = 1.0) {
  int n_time = num_meas.size();
  Vec4d rms, means, stds, p68, p95, p99;

  if (error_mat.rows() != 4) {
    std::cout << "Wrong Mat Size, Error Mat size must be (4 x timestep)" << std::endl;
    return;
  }

  // extract statistics range data
  int start_idx = (int)((1.0 - data_ratio) * n_time);
  int end_idx = n_time - 1;
  int n_range = end_idx - start_idx;

  VecXd num_meas_range(n_range);
  MatXd error_mat_range(4, n_range);

  num_meas_range = num_meas.segment(start_idx, n_range);
  error_mat_range = error_mat.block(0, start_idx, 4, n_range);

  // compute statistics ----------------------------------------------------
  // rms
  for (int i = 0; i < 4; i++) {
    rms(i) = RootMeanSquareD(error_mat_range.row(i));
    means(i) = error_mat_range.row(i).mean();
    stds(i) = StdD(error_mat_range.row(i));
    p68(i) = PercentileD(error_mat_range.row(i), 0.68);
    p95(i) = PercentileD(error_mat_range.row(i), 0.95);
    p99(i) = PercentileD(error_mat_range.row(i), 0.99);
  }

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "< Simulation Statistics (Last " << data_ratio * 100 << "%)>" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Statistics  | Position [m]  | Velocity [mm/s] | Clock Bias [ns] "
               "| Clk Drift [ns/s] "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------------"
            << std::endl;

  std::cout.precision(5);
  std::cout << "RMS         | " << std::left << std::setw(16) << rms(0) << "  " << std::left
            << std::setw(16) << rms(1) << "   " << std::left << std::setw(16) << rms(2) << std::left
            << std::setw(16) << rms(3) << std::endl;
  std::cout << "Mean+-Std   | " << std::left << means(0) << "+-" << std::left << stds(0) << "  "
            << std::left << means(1) << "+-" << std::left << stds(1) << "   " << std::left
            << means(2) << "+-" << std::left << stds(2) << "   " << std::left << means(3) << "+-"
            << std::left << stds(3) << std::endl;
  std::cout << "68%         | " << std::left << std::setw(16) << p68(0) << "  " << std::left
            << std::setw(16) << p68(1) << "   " << std::left << std::setw(16) << p68(2) << std::left
            << std::setw(16) << p68(3) << std::endl;
  std::cout << "95%         | " << std::left << std::setw(16) << p95(0) << "  " << std::left
            << std::setw(16) << p95(1) << "   " << std::left << std::setw(16) << p95(2) << std::left
            << std::setw(16) << p95(3) << std::endl;
  std::cout << "99%         | " << std::left << std::setw(16) << p99(0) << "  " << std::left
            << std::setw(16) << p99(1) << "   " << std::left << std::setw(16) << p99(2) << std::left
            << std::setw(16) << p99(3) << std::endl;
  std::cout << "  " << std::endl;
}

void Plot3DTrajectory(const std::shared_ptr<DataHistory> data_history,
                      std::string state_type = "true") {
  using namespace matplot;

  // Plot trajectory
  auto fig = figure(true);
  auto ax = fig->current_axes();
  ax->hold(on);
  ax->grid(on);
  ax->xlabel("x [km]");
  ax->ylabel("y [km]");
  ax->zlabel("z [km]");

  std::vector<Timestamped<VecXd>> data;
  if (state_type == "true") {
    ax->title("True Trajectory");
    data = data_history->GetData("rv");
  } else if (state_type == "est") {
    ax->title("Estimated Trajectory");
    data = data_history->GetData("rv_est");
  } else {
    std::cout << "State type must be either true or est." << std::endl;
    return;
  }

  // get x, y, z out of rv
  std::vector<double> x, y, z;
  for (auto& d : data) {
    x.push_back(d.GetData()(0));
    y.push_back(d.GetData()(1));
    z.push_back(d.GetData()(2));
  }
  ax->plot3(x, y, z, "b");

  // show figure
  show();
};

void PlotState(const std::shared_ptr<DataHistory> data_history, std::string state_type = "true") {
  using namespace matplot;

  if (state_type != "true" && state_type != "est") {
    std::cout << "State type must be either true or est." << std::endl;
    return;
  }

  std::vector<Timestamped<VecXd>> data_rv;
  std::vector<Timestamped<VecXd>> data_clk;

  if (state_type == "true") {
    data_rv = data_history->GetData("rv");
    data_clk = data_history->GetData("clk");
  } else {
    data_rv = data_history->GetData("rv_est");
    data_clk = data_history->GetData("clk_est");
  }

  // Plot the time history of each state, position in subplot 1, velocity in 2,
  // clock in 3
  auto fig = figure(true);
  auto ax1 = subplot(3, 1, 0);
  auto ax2 = subplot(3, 1, 1);
  auto ax3 = subplot(3, 1, 2);

  ax1->hold(on);
  ax1->grid(on);
  ax1->xlabel("Time [s]");
  ax1->ylabel("Position [km]");
  if (state_type == "true") {
    ax1->title("True Position");
  } else {
    ax1->title("Estimated Position");
  }

  ax2->hold(on);
  ax2->grid(on);
  ax2->xlabel("Time [s]");
  ax2->ylabel("Velocity [km/s]");
  if (state_type == "true") {
    ax2->title("True Velocity");
  } else {
    ax2->title("Estimated Velocity");
  }

  ax3->hold(on);
  ax3->grid(on);
  ax3->xlabel("Time [s]");
  ax3->ylabel("Clock [s]");
  if (state_type == "true") {
    ax3->title("True Clock");
  } else {
    ax3->title("Estimated Clock");
  }

  // get x, y, z out of rv
  std::vector<double> t, x, y, z, vx, vy, vz, clk_bias, clk_drift;
  for (auto& d : data_rv) {
    t.push_back(d.GetTimestamp());
    x.push_back(d.GetData()(0));
    y.push_back(d.GetData()(1));
    z.push_back(d.GetData()(2));
    vx.push_back(d.GetData()(3));
    vy.push_back(d.GetData()(4));
    vz.push_back(d.GetData()(5));
  }
  for (auto& d : data_clk) {
    clk_bias.push_back(d.GetData()(0));
    clk_drift.push_back(d.GetData()(1));
  }

  ax1->plot(t, x, "b");
  ax1->plot(t, y, "r");
  ax1->plot(t, z, "g");

  ax2->plot(t, vx, "b");
  ax2->plot(t, vy, "r");
  ax2->plot(t, vz, "g");

  ax3->plot(t, clk_bias, "b");
  ax3->plot(t, clk_drift, "r");

  ax1->legend({"x [km]", "y [km]", "z [km]"});
  ax2->legend({"vx [km/s]", "vy [km/s]", "vz [km/s]"});
  ax3->legend({"clk_bias [s]", "clk_drift [s/s]"});
  show();
};

int main() {
  /**********************************************
   * Simulation Parameters
   *********************************************/
  // Get initial epoch
  auto gps_const = GnssConstellation();
  double epoch0 = gps_const.GetEpoch();

  // Time
  double t0 = epoch0;
  double dt = 1.0;  // Integration time step [s]
  double Dt = 5.0;  // Propagation time step [s]  (= Measurement time step)
  double print_every = 600;
  double save_every = Dt;

  // Simulation seed
  int seed = 1;
  std::srand(seed);

  // Initial State
  Real a = 6541.4;
  Real e = 0.6;
  Real i = 65.5 * RAD;
  Real Omega = 0.0 * RAD;
  Real w = 90.0 * RAD;
  Real M = 0.0 * RAD;
  Real clk_bias = 0.0;
  Real clk_drift = 0.1;

  // Set simulation to 1 orbit
  int n_orbit = 2;  // number of orbits to simulate
  Real period = 2.0 * M_PI * sqrt(pow(a, 3) / GM_MOON);
  double tf = t0 + n_orbit * period.val();
  int time_step_num = int((tf - t0) / Dt) + 1;
  tf = t0 + (time_step_num - 1) * Dt;

  // Dynamics Model   Todo: Refine this to a more high fidelity model
  int moon_sph_true = 8;  // moon spherical harmonics order in true dynamics
  int moon_sph_est = 5;   // moon spherical harmonics order in filter dynamics
  bool add_earth = true;  // add earth to true and filter dynamics

  // Onboard Clock Model
  ClockModel cmodel = ClockModel::kMiniRafs;

  // measurements
  bool use_range = true;       // use GPS pseudorange measurement
  bool use_range_rate = true;  // use GPS pseudorange-rate measurement

  // Estimation
  int state_size = 8;           // Pos(3), vel(3), bias, drift [km, km/s, s, s/s]
  double pos_err = 1.0;         // Initial Position error [km]
  double vel_err = 1e-3;        // Initial Velocity error [km/s]
  double clk_bias_err = 1e-6;   // Initial Clock bias error [s]
  double clk_drift_err = 1e-9;  // Initial Clock drift error [s/s]
  double sigma_acc = 1e-12;     // Process noise Acceleration [km/s^2]  <-- tune
                                // this for optimal performance!

  // Debug mode
  bool plot_results = false;
  bool debug_jacobian = false;
  bool print_debug = false;
  bool debug_ekf = false;
  bool debug_ekf_error_only = true;  // Print only error for EKF debugging
  bool no_meas = false;              // set to true to turn off measurements

  /**********************************************
   * Setup
   * ********************************************/

  // measurement type vector
  std::vector<GnssMeasurementType> meas_types;
  if (use_range) {
    meas_types.push_back(GnssMeasurementType::PR);
  }
  if (use_range_rate) {
    meas_types.push_back(GnssMeasurementType::PRR);
  }
  int meas_type_num = meas_types.size();

  // Orbit Dynamics
  IntegratorParams iparams;
  iparams.abstol = 1e-12;
  iparams.reltol = 1e-12;

  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(
      GM_EARTH);  // use 2d earth dynamics to propagate GPS constellation
  auto dyn_est = MakePtr<NBodyDynamics<Real>>(IntegratorType::RKF45);     // Filter Dynamics
  auto dyn_true = MakePtr<NBodyDynamics<double>>(IntegratorType::RKF45);  // true dynamics

  dyn_true->SetIntegratorParams(iparams);
  dyn_est->SetIntegratorParams(iparams);

  auto earth_true = BodyT<double>::Earth();
  auto earth_est = BodyT<Real>::Earth();
  auto moon_true = BodyT<double>::Moon(moon_sph_true, moon_sph_true);
  auto moon_est = BodyT<Real>::Moon(moon_sph_est, moon_sph_est);

  dyn_true->SetFrame(Frame::MOON_CI);
  dyn_est->SetFrame(Frame::MOON_CI);
  dyn_true->AddBody(moon_true);
  dyn_est->AddBody(moon_est);

  if (add_earth) {
    dyn_true->AddBody(earth_true);
    dyn_est->AddBody(earth_est);
  }

  // clock dynamics
  auto dyn_clk_true = ClockDynamics(cmodel);
  auto dyn_clk_est = ClockDynamics(cmodel);
  dyn_clk_true.SetNoise(true);
  dyn_clk_est.SetNoise(false);

  // GPS constellation
  auto channel = std::make_shared<GnssChannel>();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");  // example gps file

  // Print Time
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  // Set dynamics integration time
  dyn_earth_tb->SetTimeStep(dt);
  dyn_est->SetTimeStep(dt);
  dyn_true->SetTimeStep(dt);

  // Moon spacecraft
  ClassicalOE coe_moon({a, e, i, Omega, w, M}, Frame::MOON_CI);
  auto cart_state_moon = std::make_shared<CartesianOrbitState>(Classical2Cart(coe_moon, GM_MOON));

  Vec2 clock_vec{clk_bias, clk_drift};  // [s, s/s]
  ClockState clock_state(clock_vec);

  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GnssReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_true);
  moon_sat->SetClock(clock_state);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(NaifId::MOON);
  moon_sat->SetClockDynamics(dyn_clk_true);

  receiver->SetAgent(moon_sat);
  receiver->SetReceiverAttitudeMode("PZ_EarthPoint");
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Initial covariance
  MatXd P0 = ConstructInitCovariance(pos_err, vel_err, clk_bias_err, clk_drift_err);

  // Joint state and dynamics
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(cart_state_moon.get(), dyn_est.get());
  joint_state.PushBackStateAndDynamics(&clock_state, &dyn_clk_est);

  FilterDynamicsFunction joint_dynamics = joint_state.GetFilterDynamicsFunction();

  /*********************************************
   * Define Measurement function
   * *******************************************/
  FilterMeasurementFunction meas_func_pos_clk
      = [moon_sat, receiver, state_size, no_meas, meas_types, debug_jacobian](
            const VecX x, MatXd& H, MatXd& R) -> VecX {
    if (no_meas) {
      return VecXd::Zero(0);
    }

    // Measurements
    std::string freq = "L1";
    double epoch = moon_sat->GetEpoch().val();
    auto measall = receiver->GetMeasurement(epoch);  // measurements of all frequencies
    auto meas = measall.ExtractSignal("L1");         // measurements of L1
    int sat_num = meas.GetTrackedSatelliteNum();     // number of tracked GPS satellites
    int mtot = sat_num * meas_types.size();          // total number of measurements

    // Predict measurements
    H = MatXd::Zero(mtot, x.size());
    R = MatXd::Zero(mtot, mtot);
    VecX x_N = VecX::Zero(mtot);  // a dummy variable for carrier phase
    Frame frame_in = Frame::MOON_CI;

    VecX z = meas.GetPredictedGnssMeasurement(epoch, x.head(6), x.tail(2), x_N, H, meas_types,
                                              frame_in);  // Jacobian with autodiff

    if (debug_jacobian) {
      // Compute Numerical Jacobian
      MatXd H_num = MatXd::Zero(mtot, x.size());
      MatXd H_dum = MatXd::Zero(mtot, x.size());

      for (int i = 0; i < x.size(); i++) {
        Real eps = x(i) * 1e-6;
        VecX x_p = x;
        VecX x_m = x;
        x_p(i) += eps;
        x_m(i) -= eps;
        VecX z_p = meas.GetPredictedGnssMeasurement(epoch, x_p.head(6), x_p.tail(2), x_N, H_dum,
                                                    meas_types,
                                                    frame_in);  // Jacobian with autodiff
        VecX z_m = meas.GetPredictedGnssMeasurement(epoch, x_m.head(6), x_m.tail(2), x_N, H_dum,
                                                    meas_types,
                                                    frame_in);  // Jacobian with autodiff
        H_num.col(i) = ((z_p - z_m) / (2 * eps)).cast<double>();
      }

      std::cout << " " << std::endl;
      std::cout << "AutoDiff Jacobian: " << std::endl << H << std::endl;
      std::cout << "Numerical Jacobian: " << std::endl << H_num << std::endl;
      std::cout << "Jacobian Error: " << (H - H_num).norm() << std::endl;
      std::cout << " " << std::endl;
    }

    // Get the Measurement Noise
    VecXd noise_std_vec = meas.GetGnssNoiseStdVec(meas_types);
    R.diagonal().array() = noise_std_vec.array().square();

    return z;
  };

  /*********************************************
   * Define Process Noise function
   * *******************************************/
  FilterProcessNoiseFunction proc_noise_func
      = [cmodel, state_size, sigma_acc](const VecX x, Real t_curr, Real t_end) -> MatXd {
    int clock_index = 6;
    double dt = (t_end - t_curr).val();

    MatXd Q = MatXd::Zero(state_size, state_size);

    Mat6d Q_rv = Mat6d::Zero();
    for (int i = 0; i < 3; i++) {
      Q_rv(i, i) = pow(dt, 3) / 3.0 * pow(sigma_acc, 2);
      Q_rv(i + 3, i + 3) = dt * pow(sigma_acc, 2);
      Q_rv(i, i + 3) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
      Q_rv(i + 3, i) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
    }

    Mat2d Q_clk = ClockDynamics::TwoStateNoise(cmodel, dt).cast<double>();

    Q.block(0, 0, 6, 6) = Q_rv;
    Q.block(6, 6, 2, 2) = Q_clk;

    return Q;
  };

  /*************************************
   * EKF Setup
   * ***********************************/
  EKF ekf;
  ekf.SetDynamicsFunction(joint_dynamics);
  ekf.SetMeasurementFunction(meas_func_pos_clk);
  ekf.SetProcessNoiseFunction(proc_noise_func);
  std::cout << "Initialized EKF" << std::endl;

  // Storage
  MatXd error_mat(4, time_step_num);
  VecXd num_meas(time_step_num);

  // Initilization
  VecX x_est = SampleMVN(joint_state.GetJointStateValue(), P0, 1, seed);
  ekf.Initialize(x_est, P0);
  VecXd est_err = ComputeEstimationErrors(moon_sat, &ekf);
  error_mat.col(0) = est_err;

  // Print State
  if (print_debug) {
    VecX x_init_true = moon_sat->GetStateVec();

    std::cout << "Initial true state: " << moon_sat->GetStateVec().transpose() << std::endl;
    std::cout << "Initial estimated state: " << x_est.transpose() << std::endl;
    VecX err_std = x_est - x_init_true;
    for (int i = 0; i < 8; i++) {
      err_std(i) = err_std(i) / sqrt(P0(i, i));
    }
    std::cout << "Initial error / std: " << err_std << std::endl;
  }

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = std::filesystem::current_path() / "output" / "ExampleEKF";
  FileWriter writer(output_path, true);

  /***********************************************
   * Main loop
   **********************************************/
  Real t = t0;

  double epoch = epoch0;
  PrintProgressHeader();

  int time_index = 0;
  // tf = 50 * Dt;

  // Compute Estimation
  est_err = ComputeEstimationErrors(moon_sat, &ekf);  // pos, vel, clkb, clkd error
  PrintProgress(t.val(), est_err(0), est_err(1), est_err(2));
  for (t = t0; t < tf; t += Dt) {
    time_index += 1;
    epoch += Dt;  // first propagate to the next epoch

    // Propagate True State
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    // Get True Measurement
    VecX z_true;
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int num_sat = meas.GetTrackedSatelliteNum();
    num_meas(time_index) = num_sat;

    if (!no_meas) {
      bool with_noise = true;
      z_true = meas.GetGnssMeasurement(meas_types, with_noise, seed);
    } else {
      z_true = VecXd::Zero(0);
    }

    // Update EKF
    ekf.Predict(t + Dt);
    // print Phi
    // std::cout << "Phi:" << std::endl << ekf.Phi_ << std::endl;

    ekf.Update(z_true, debug_ekf);

    // Add Data
    if (!no_meas) {
      AddStateEstimationData(data_history, moon_sat, &ekf, &gps_const, &meas, t.val(), epoch);
    }

    // Compute Estimation
    est_err = ComputeEstimationErrors(moon_sat, &ekf);  // pos, vel, clkb, clkd error
    error_mat.col(time_index) = est_err;

    // Print progress
    if (fmod(t.val(), print_every) < 1e-3) {
      PrintProgress(t.val(), est_err(0), est_err(1), est_err(2));
      PrintEKFDebugInfo(time_index, moon_sat, &ekf, true);
    }

    // print measurement residuals
    if ((print_debug) && (!no_meas)) {
      PrintEKFDebugInfo(time_index, moon_sat, &ekf, debug_ekf_error_only);
    }
  }

  // Print Statistics
  if (!no_meas) {
    PrintEstimationStatistics(num_meas, error_mat, 0.3);  // use last 30%

    // Write data ------------------------------------------
    std::cout << "Simulation finished, saving data..." << std::endl;
    writer.WriteData(*data_history);
    std::cout << "Data saved to " << output_path << std::endl;

    // Plots -----------------------------------------------
    if (plot_results) {
      Plot3DTrajectory(data_history, "true");
      Plot3DTrajectory(data_history, "est");

      // Plot estimate for all states
      PlotState(data_history, "true");
      PlotState(data_history, "est");
    }
  }
}
