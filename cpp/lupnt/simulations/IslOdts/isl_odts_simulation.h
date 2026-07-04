#pragma once

#include <string>
#include <vector>

#include "lupnt/core/simulation.h"

namespace lupnt {

  /// @brief Initial Cartesian state and clock state for one satellite in an
  /// `IslOdtsSimulation`, defined at `IslOdtsConfig::start_epoch_utc` in the
  /// Moon-centered inertial frame (`Frame::MOON_CI`).
  struct IslOdtsSatelliteConfig {
    std::string name = "SAT";

    /// Initial Cartesian position [m] and velocity [m/s] in `Frame::MOON_CI`, e.g.
    /// taken directly from a published reference-constellation ephemeris (ICRF is
    /// treated as equivalent to `Frame::MOON_CI` for this purpose).
    Vec3d r0_m = Vec3d::Zero();
    Vec3d v0_mps = Vec3d::Zero();

    double clock_bias_s = 0.0;
    double clock_drift_sps = 0.0;
  };

  /// @brief Simplified inter-satellite-link (ISL) crosslink budget, used to report a
  /// representative carrier-to-noise density ratio alongside the range/Doppler
  /// measurements. Not used by the estimator; provided so link performance can be
  /// assessed without a full antenna-gain-pattern/transmitter model (see
  /// `lupnt::GnssMeasurement`/`lupnt::Antenna` for that level of fidelity).
  struct IslLinkBudgetConfig {
    bool enabled = true;
    double tx_power_dbw = 5.0;             // ~3 W crosslink transmit power
    double tx_gain_dbi = 20.0;             // representative crosslink antenna gain
    double rx_gain_dbi = 20.0;
    double frequency_hz = 25.5e9;          // representative Ka-band crosslink carrier
    double system_noise_temp_k = 500.0;
  };

  struct IslOdtsConfig {
    int seed = 42;
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double duration_s = 6.0 * 3600.0;
    double dt_s = 60.0;
    double integration_step_s = 60.0;

    /// Constellation of N >= 2 satellites. `satellites[0]` is the "own"/hub
    /// satellite whose onboard `SchmidtEKF` is run; `satellites[1..N-1]` are the
    /// satellites it maintains a two-way crosslink to (N-1 links total), each
    /// carried as one Schmidt consider-state block, in the given order.
    std::vector<IslOdtsSatelliteConfig> satellites;

    IslLinkBudgetConfig link_budget;

    // Force model
    int moon_gravity_degree_truth = 20;
    int moon_gravity_order_truth = 20;
    int moon_gravity_degree_filter = 8;
    int moon_gravity_order_filter = 8;
    bool include_earth = true;
    bool include_sun = true;
    bool use_relativity = true;

    // Measurement noise (two-way range / Doppler collapsed to the geometric
    // relative range and range-rate, see `isl_odts_simulation.cc`), applied to
    // every crosslink identically.
    double range_sigma_m = 1.0;
    double range_rate_sigma_mps = 1.0e-3;

    // Onboard Schmidt-EKF tuning. The hub filter estimates its own 8-state
    // [r, v, clock_bias, clock_drift] and considers each linked satellite's
    // 8-state broadcast/prior state without correcting it.
    double initial_position_sigma_m = 200.0;
    double initial_velocity_sigma_mps = 0.1;
    double initial_clock_bias_sigma_s = 1.0e-6;
    double initial_clock_drift_sigma_sps = 1.0e-9;

    // A satellite's knowledge of a neighbor's state -- e.g. from a periodically
    // refreshed broadcast/prior ephemeris -- is usually somewhat better than its
    // own just-initialized uncertainty but still not directly correctable onboard;
    // these sigmas seed the (fixed-mean) consider-state covariance block, applied
    // identically to every linked satellite.
    double consider_position_sigma_m = 100.0;
    double consider_velocity_sigma_mps = 0.05;
    double consider_clock_bias_sigma_s = 1.0e-7;
    double consider_clock_drift_sigma_sps = 1.0e-10;

    double process_accel_sigma_mps2 = 1.0e-8;
  };

  /// @brief Time series results of an `IslOdtsSimulation` run, in SI units,
  /// `Frame::MOON_CI`. Row `k` of every matrix corresponds to `t_s(k)`; column `i`
  /// of every per-link matrix corresponds to the crosslink to
  /// `satellite_names[i + 1]` (i.e. `IslOdtsConfig::satellites[i + 1]`).
  struct IslOdtsResults {
    VecXd t_s;  // Elapsed time since start_epoch_utc [s], size [N]

    std::vector<std::string> satellite_names;  // size [n_sat], [0] = own/hub satellite

    // Truth trajectories, one per satellite (same order as satellite_names), columns
    // [r_x,r_y,r_z,v_x,v_y,v_z,clock_bias_s,clock_drift_sps]
    std::vector<MatXd> truth_states;  // size [n_sat], each [N x 8]

    // Onboard Schmidt-EKF estimate/covariance-diagonal, columns
    // [own(8), consider_1(8), ..., consider_{n_sat-1}(8)]
    MatXd est;       // [N x 8*n_sat]
    MatXd cov_diag;  // [N x 8*n_sat]

    // Two-way range/Doppler (relative range-rate) measurements, one column per
    // crosslink (n_links = n_sat - 1)
    MatXd range_true_m;         // [N x n_links]
    MatXd range_rate_true_mps;  // [N x n_links]
    MatXd range_obs_m;          // [N x n_links]
    MatXd range_rate_obs_mps;   // [N x n_links]

    // Pre-fit measurement residuals (observed minus the hub filter's predicted
    // value, evaluated just before that epoch's Update() corrects the state) -- the
    // standard OD diagnostic for how well the filter is tracking each crosslink,
    // independent of the (possibly weakly observable) absolute position error.
    MatXd range_resid_m;         // [N x n_links]
    MatXd range_rate_resid_mps;  // [N x n_links]

    // Simplified crosslink carrier-to-noise density ratio [dB-Hz] (NaN if disabled)
    MatXd cn0_dbhz;  // [N x n_links]
  };

  /// @brief Onboard inter-satellite-link (ISL) orbit determination and timing system
  /// (ODTS) simulation for a constellation of N satellites (e.g. Elliptical Lunar
  /// Frozen Orbit satellites from a lunar relay/navigation constellation), where one
  /// "hub" satellite maintains a two-way crosslink to each of the other N-1.
  ///
  /// Propagates a truth trajectory for every satellite, simulates two-way range and
  /// Doppler (range-rate) crosslink measurements between the hub and each other
  /// satellite, and runs one Schmidt Extended Kalman Filter (`SchmidtEKF`,
  /// `lupnt/filters/schmidt_ekf.h`) onboard the hub. The filter estimates the hub's
  /// own `[r, v, clock_bias, clock_drift]` state from the crosslink measurements
  /// while carrying each other satellite's state as a separate Schmidt "consider"
  /// state: its uncertainty (and correlation with the estimated state) is propagated
  /// and used in the measurement update, but it is never corrected, mirroring a hub
  /// that only has a broadcast/prior estimate of its neighbors' states.
  class IslOdtsSimulation : public Simulation {
  public:
    explicit IslOdtsSimulation(IslOdtsConfig config);

    void Setup() override;
    void Run() override;

    const IslOdtsConfig& GetConfig() const { return config_; }
    const IslOdtsResults& GetResults() const { return results_; }

  private:
    IslOdtsConfig config_;
    IslOdtsResults results_;
    bool setup_complete_ = false;
  };

}  // namespace lupnt
