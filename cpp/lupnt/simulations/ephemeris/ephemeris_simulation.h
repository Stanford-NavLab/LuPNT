#pragma once

#include <string>
#include <vector>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/simulations/simulation.h"

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
    /// Number of Fourier terms (harmonics of the argument of latitude) added to
    /// the ephemeris Chebyshev residual model; 0 (default) is a pure Chebyshev
    /// model. See `EphemerisFitOptions::num_fourier_terms`.
    int cartesian_num_fourier_terms = 0;
    int almanac_poly_order = 1;
    /// Number of Fourier harmonics per element in the almanac model (see
    /// `AlmanacFitOptions::num_fourier_terms`).
    int almanac_num_fourier_terms = 1;

    /// Frame the ephemeris/almanac are fit and evaluated in. Defaults to
    /// `propagate_frame` (inertial). Set to `Frame::MOON_PA` to fit and broadcast
    /// in the rotating Moon-fixed principal-axis frame: the truth trajectory is
    /// converted from `propagate_frame` to this frame before fitting (per Iiyama &
    /// Gao). When equal to `propagate_frame` no conversion is applied.
    Frame output_frame = Frame::MOON_CI;

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

  /// @brief Full result payload of an ephemeris datasize/accuracy sweep, exposed by
  /// `EphemerisApp::GetResults()` (identical layout to the former
  /// `EphemerisSimulation`'s individual accessors). All series are in SI units, in
  /// `EphemerisSimulationConfig::output_frame`.
  struct EphemerisResults {
    VecXd t_truth_s;  // Elapsed time since start_epoch_utc [s], size [N]
    MatXd rv_truth;   // Truth Cartesian states [N x 6] in output_frame

    /// One entry per `EphemerisSimulationConfig::fit_window_minutes` value.
    std::vector<EphemerisWindowResult> cartesian_results;
    std::vector<EphemerisWindowResult> almanac_results;
  };

  // The ephemeris/almanac datasize-accuracy study is driven by `EphemerisApp`
  // (applications/ephemeris/ephemeris_app.h) on a thin `EphemerisManager` agent, via
  // `pnt.Simulation(...)`. The `EphemerisSimulationConfig`/`EphemerisResults` structs above
  // are the shared config/result payloads consumed by that app: for each
  // `fit_window_minutes` entry, `num_windows` windows are sampled across a numerically
  // propagated lunar-orbit truth trajectory; each ephemeris/almanac model is fit and
  // evaluated (RMS/95th-percentile RTN position and velocity error) on every window, and the
  // minimum per-parameter broadcast resolution (bit count) needed to keep the position
  // accuracy within `datasize_precision_m` is found by bisection -- mirroring the
  // accuracy/datasize trade-off studies used to size GNSS broadcast ephemeris/almanac
  // message formats.

}  // namespace lupnt
