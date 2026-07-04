#pragma once

#include <string>
#include <vector>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/simulation.h"

namespace lupnt {

  /// @brief Initial classical orbital elements for the truth trajectory sampled by
  /// `EphemerisSimulation` (e.g. a lunar Elliptical Frozen Orbit).
  struct EphemerisOrbitConfig {
    double a_m = 6541.4e3;
    double ecc = 0.6;
    double inc_rad = 56.2 * RAD;
    double raan_rad = 0.0;
    double argp_rad = 90.0 * RAD;
    double m0_rad = 0.0;

    /// Frame the elements above are defined in (e.g. `Frame::MOON_OP` for a
    /// frozen-orbit definition); converted to `EphemerisSimulationConfig::
    /// propagate_frame` at `start_epoch_utc`.
    Frame coe_frame = Frame::MOON_OP;
  };

  struct EphemerisSimulationConfig {
    std::string start_epoch_utc = "2027-01-01T00:00:00";
    EphemerisOrbitConfig orbit;

    /// Frame the truth trajectory is propagated and ephemerides/almanacs are fit
    /// in (e.g. `Frame::MOON_CI`).
    Frame propagate_frame = Frame::MOON_CI;

    double duration_days = 3.0;
    double sample_dt_s = 60.0;

    // Force model for the numerically propagated truth trajectory.
    int moon_gravity_degree = 20;
    int moon_gravity_order = 20;
    bool include_earth = true;
    bool include_sun = true;
    bool use_relativity = false;
    double integration_step_s = 30.0;

    /// Fitting-window lengths [minutes] swept by the ephemeris/almanac
    /// accuracy-vs-datasize study.
    std::vector<double> fit_window_minutes = {60.0, 120.0, 240.0, 480.0};

    /// Number of fitting windows sampled across the truth trajectory for each
    /// `fit_window_minutes` entry; fit accuracy is averaged, and the per-parameter
    /// broadcast resolution (and hence total bit count) is maximized, across these
    /// samples.
    int num_windows = 6;

    int cartesian_poly_order = 8;
    bool cartesian_use_keplerian_baseline = true;
    int almanac_poly_order = 2;

    /// Required position accuracy [m] used to search each parameter's minimum
    /// broadcast resolution (see `EphemerisSimulation::Run`).
    double datasize_precision_m = 0.01;
  };

  /// @brief Fit-accuracy / broadcast-datasize summary for one
  /// `fit_window_minutes` entry, averaged (accuracy) or maximized (bit count)
  /// across `EphemerisSimulationConfig::num_windows` windows sampled across the
  /// orbit.
  struct EphemerisWindowResult {
    double fit_window_min = 0.0;
    int num_params = 0;
    int total_bits = 0;
    double pos_rms_m = 0.0;
    double vel_rms_mps = 0.0;
    double pos_p95_m = 0.0;
    double vel_p95_mps = 0.0;
  };

  /// @brief Studies the position/velocity fit-accuracy vs. broadcast-datasize
  /// trade-off of `CartesianEphemeris` (precise, short-validity) and `Almanac`
  /// (coarse, long-validity) orbit models -- both in `lupnt/applications/` -- fit
  /// to a numerically propagated lunar orbit truth trajectory.
  ///
  /// For each `fit_window_minutes` entry, `num_windows` windows are sampled across
  /// the truth trajectory; each ephemeris/almanac model is fit and evaluated
  /// (RMS/95th-percentile RTN position and velocity error) on every window, and
  /// the minimum per-parameter broadcast resolution (bit count) needed to keep the
  /// position accuracy within `datasize_precision_m` is found by bisection --
  /// mirroring the accuracy/datasize trade-off studies used to size GNSS broadcast
  /// ephemeris/almanac message formats.
  class EphemerisSimulation : public Simulation {
  public:
    explicit EphemerisSimulation(EphemerisSimulationConfig config);

    void Setup() override;
    void Run() override;

    const EphemerisSimulationConfig& GetConfig() const { return config_; }
    const VecXd& GetTruthTimes() const { return t_truth_s_; }
    const MatXd& GetTruthStates() const { return rv_truth_; }
    const std::vector<EphemerisWindowResult>& GetCartesianResults() const {
      return cartesian_results_;
    }
    const std::vector<EphemerisWindowResult>& GetAlmanacResults() const { return almanac_results_; }

  private:
    EphemerisSimulationConfig config_;
    bool setup_complete_ = false;

    VecXd t_truth_s_;  // Elapsed time since start_epoch_utc [s]
    MatXd rv_truth_;   // Truth Cartesian states [N x 6] in propagate_frame

    std::vector<EphemerisWindowResult> cartesian_results_;
    std::vector<EphemerisWindowResult> almanac_results_;
  };

}  // namespace lupnt
