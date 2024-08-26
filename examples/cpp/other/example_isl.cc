/**
 * @file example_isl.cc
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-08-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <lupnt/lupnt.h>

using namespace lupnt;
namespace sp = spice;

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

VecX ExtractSatState(const VecX& x, int sat_idx, int state_per_sat) {
  return x.segment(sat_idx * state_per_sat, state_per_sat);
};

class ISLTransmitter : public Transmitter {
public:
  double GetTransmitterAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) override {
    return 0.0;
  };
};

class ISLReceiver : public Receiver {
public:
  double GetReceiverAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) override {
    return 0.0;
  };
};

int main() {
  /**********************************************
   * Parameter Setting
   * ********************************************/
  // Time
  Real epoch0 = spice::String2TAI("2035/02/01 00:00:00.000 UTC");
  double t0 = epoch0.val();
  double dt = 1.0;  // Integration time step [s]
  double Dt = 5.0;  // Propagation time step [s]  (= Measurement time step)
  double print_every = 600;
  double save_every = Dt;

  // Simulation seed
  int seed = 1;
  std::srand(seed);

  int nsat = 3;  // number of satellites

  // IOAG study report
  // https://www.ioag.org/Public%20Documents/Lunar%20communications%20architecture%20study%20report%20FINAL%20v1.3.pdf
  Real M1 = 0.0;    // mean anomaly of satellite 1 (rad)
  Real M2 = 0.0;    // mean anomaly of satellite 2 (rad)
  Real M3 = 0.0;    // mean anomaly of satellite 3 (rad)
  Real a = 6142.4;  // semi-major axis (km)
  MatX oe(6, nsat);
  Vec6 sat1_oe = {a, 0.0, 0.0, 0.0, RAD * 315, M1};
  Vec6 sat2_oe = {a, 0.6, RAD * 57.7, RAD * 270, RAD * 270, M2};  // 2
  Vec6 sat3_oe = {a, 0.6, RAD * 57.7, 0.0, RAD * 90, M3};         // 3
  oe.col(0) = sat1_oe;
  oe.col(1) = sat2_oe;
  oe.col(2) = sat3_oe;

  // Set simulation to 1 orbit
  int n_orbit = 1;  // number of orbits to simulate
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
  bool use_range = true;       // use range measurements
  bool use_range_rate = true;  // use range rates measurement
  bool use_fixed_error = true;
  double range_sigma_fixed = 1e-3;
  double range_rate_sigma_fixed = 1e-6;

  // Occultations
  std::vector<NaifId> occult_bodies = {NaifId::MOON};
  VecXd occult_alt(1);
  occult_alt << 100.0;         // occultation altitude [km]
  Real hardware_delay = 1e-8;  // hardware delay [s]

  // Estimation
  int state_size = 8;           // Pos(3), vel(3), bias, drift [km, km/s, s, s/s]
  double pos_err = 1.0;         // Initial Position error [km]
  double vel_err = 1e-3;        // Initial Velocity error [km/s]
  double clk_bias_err = 1e-6;   // Initial Clock bias error [s]
  double clk_drift_err = 1e-9;  // Initial Clock drift error [s/s]
  double sigma_acc = 1e-12;     // Process noise Acceleration [km/s^2]  <-- tune
                                // this for optimal performance!

  // Transmitter
  Ptr<ISLTransmitter> transmitter = MakePtr<ISLTransmitter>();
  transmitter->P_tx = 10.0;       // Transmit power [dBW]
  transmitter->freq_tx = 26.5e9;  // Transmit frequency [Hz]

  // Receiver
  Ptr<ISLReceiver> receiver = MakePtr<ISLReceiver>();
  Modulation modulation_type = Modulation::BPSK;  // carrier type
  receiver->rx_param_.B_L_chip = 0.1;             // tracking loop noise bandwidth [Hz]
  receiver->rx_param_.Tc = 1 / 2.068e6;           // chip duration
  receiver->rx_param_.B_L_carrier = 0.1;          // carrier loop noise bandwidth [Hz]
  receiver->rx_param_.m_R = 0.0;                  // modulation index
  receiver->rx_param_.T_I_doppler = 10.0;         // Doppler integration time [s]
  receiver->rx_param_.T_I_range = 0.5;            // range integration time [s] (for open loop)
  receiver->rx_param_.pn_ranging_code = "T4B";    //  "T2B", "T4B"
  receiver->rx_param_.SER_threshold = 0.1;        // Symbol error rate threshold
  receiver->rx_param_.BTs = 0.5;  // For GMSK modulation (B: 3dB point of gaussian filter, )
  receiver->rx_param_.coding_rate = 1 / 2;  // Coding rate

  /**********************************************
   * Setup
   * ********************************************/

  // measurement type vector
  std::vector<LinkMeasurementType> meas_types;
  if (use_range) {
    meas_types.push_back(LinkMeasurementType::Range);
  }
  if (use_range_rate) {
    meas_types.push_back(LinkMeasurementType::RangeRate);
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

  auto moon_true = BodyT<double>::Moon(moon_sph_true, moon_sph_true);
  auto moon_est = BodyT<Real>::Moon(moon_sph_est, moon_sph_est);

  dyn_true->SetFrame(Frame::MOON_CI);
  dyn_est->SetFrame(Frame::MOON_CI);
  dyn_true->AddBody(moon_true);
  dyn_est->AddBody(moon_est);

  if (add_earth) {
    dyn_true->AddBody(BodyT<double>::Earth());
    dyn_est->AddBody(BodyT<Real>::Earth());
  }

  // clock dynamics
  auto dyn_clk_true = ClockDynamics(cmodel);
  auto dyn_clk_est = ClockDynamics(cmodel);
  dyn_clk_true.SetNoise(true);
  dyn_clk_est.SetNoise(false);

  // Print Time
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  // Set dynamics integration time
  dyn_earth_tb->SetTimeStep(dt);
  dyn_est->SetTimeStep(dt);
  dyn_true->SetTimeStep(dt);

  // Moon spacecraft
  std::vector<Ptr<Spacecraft>> moon_sats;
  JointState joint_state;

  for (int i = 0; i < nsat; i++) {
    // orbit
    ClassicalOE coe_moon(oe.col(i), Frame::MOON_CI);
    Ptr<CartesianOrbitState> cart_state
        = MakePtr<CartesianOrbitState>(Classical2Cart(coe_moon, GM_MOON));
    Ptr<Spacecraft> moon_sat = MakePtr<Spacecraft>();

    // clock
    Real clk_bias = Real(SampleRandNormal(0.0, clk_bias_err, seed));
    Real clk_drift = Real(SampleRandNormal(0.0, clk_drift_err, seed));
    Vec2 clock_vec{clk_bias, clk_drift};  // [s, s/s]
    ClockState clock_state(clock_vec);

    auto transponder = MakePtr<Transponder>(transmitter, receiver);

    moon_sat->AddDevice(transponder);
    moon_sat->SetDynamics(dyn_true);
    moon_sat->SetClock(clock_state);
    moon_sat->SetOrbitState(cart_state);
    moon_sat->SetEpoch(epoch0);
    moon_sat->SetBodyId(NaifId::MOON);
    moon_sat->SetClockDynamics(dyn_clk_true);

    joint_state.PushBackStateAndDynamics(cart_state.get(), dyn_est.get());
    joint_state.PushBackStateAndDynamics(&clock_state, &dyn_clk_est);
  }

  // Initial covariance
  MatXd P0 = ConstructInitCovariance(pos_err, vel_err, clk_bias_err, clk_drift_err);

  FilterDynamicsFunction joint_dynamics = joint_state.GetFilterDynamicsFunction();

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

  /*********************************************
   * Define Measurement function
   * *******************************************/
  std::vector<std::pair<int, int>> sat_pairs_idx;
  bool no_meas = false;
  double epoch_rx = 0;

  LinkMeasurement link_meas = LinkMeasurement(occult_bodies, occult_alt, hardware_delay);
  if (use_fixed_error) {
    link_meas.UseFixedError();
    link_meas.SetFixedRangeError(range_sigma_fixed);
    link_meas.SetFixedRangeRateError(range_rate_sigma_fixed);
  }

  FilterMeasurementFunction meas_func_pos_clk
      = [link_meas, epoch_rx, sat_pairs_idx, state_size, no_meas, meas_types](
            const VecX x, MatXd& H, MatXd& R) -> VecX {
    if (no_meas) {
      return VecX::Zero(0);
    }

    // Total Number of Measurements
    int mtot = meas_types.size() * sat_pairs_idx.size();
    H.resize(mtot, x.size());
    R.resize(mtot, mtot);
    VecX z = VecX::Zero(mtot);

    // Iterate over all ISL pairs
    for (int i = 0; i < sat_pairs_idx.size(); i++) {
      int sat_target_idx = sat_pairs_idx[i].first;
      int sat_rx_idx = sat_pairs_idx[i].second;

      // Measurements
      VecX sat_target = ExtractSatState(x, sat_target_idx, state_size);
      VecX sat_rx = ExtractSatState(x, sat_rx_idx, state_size);

      int mtot = meas_types.size();

      // Predict Measurements
      // VecX z = link_meas.GetPredictedGnssMeasurement(
      //     epoch_, x.head(6), x.tail(2), x_N, H, meas_types,
      //     frame_in);  // Jacobian with autodif

      // // Get the Measurement Noise
      // VecXd noise_std_vec = meas.GetGnssNoiseStdVec(meas_types);
      // R.diagonal().array() = noise_std_vec.array().square();
    }
    return z;
  };
}
