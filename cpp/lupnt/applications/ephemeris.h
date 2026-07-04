#pragma once

#include <string>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Position/velocity fit-error summary statistics for an ephemeris/almanac
  /// evaluation, decomposed into Radial/Transverse/Normal (RTN) components at the
  /// reference trajectory, plus the combined 3-D norm (index 3).
  struct EphemerisFitErrorStats {
    Vec4d rms_pos_m = Vec4d::Zero();    // [R, T, N, 3D] position RMS error [m]
    Vec4d rms_vel_mps = Vec4d::Zero();  // [R, T, N, 3D] velocity RMS error [m/s]
    Vec4d p95_pos_m = Vec4d::Zero();    // [R, T, N, 3D] 95th-percentile |position error| [m]
    Vec4d p95_vel_mps = Vec4d::Zero();  // [R, T, N, 3D] 95th-percentile |velocity error| [m/s]
  };

  /// @brief Configuration for `CartesianEphemeris::Fit`/`Eval`.
  struct EphemerisFitOptions {
    /// Order of the Chebyshev polynomial used to represent the position residual
    /// (and, through its analytic time derivative, the velocity residual) left
    /// over after subtracting the two-body Keplerian baseline orbit.
    int poly_order = 8;

    /// If true, subtract an osculating two-body Kepler orbit (fit at the fitting
    /// window's midpoint epoch) before fitting the Chebyshev polynomial to the
    /// residual -- this is what lets a modest polynomial order stay accurate over
    /// a multi-hour window. If false, position/velocity are represented by the
    /// polynomial alone (useful as a baseline comparison, but needs a much higher
    /// order for the same accuracy).
    bool use_keplerian_baseline = true;

    /// Gravitational parameter of the central body [m^3/s^2], consistent with the
    /// position/velocity units passed to Fit/Eval.
    double gm = GM_MOON;
  };

  /// @brief Piecewise ephemeris model for one fitting window: an osculating
  /// two-body Kepler orbit (fit at the window's midpoint) plus a
  /// Chebyshev-polynomial correction on the Cartesian position residual, with
  /// velocity obtained from the polynomial's analytic time derivative. This
  /// mirrors the "orbital elements + polynomial correction" representation used by
  /// GNSS broadcast ephemerides (e.g. GPS LNAV), extended here to arbitrary
  /// orbits/gravity fields via the polynomial correction term.
  ///
  /// Used directly (via the Python bindings) and by `EphemerisSimulation`
  /// (`lupnt/simulations/Ephemeris/ephemeris_simulation.h`) to study the
  /// position/velocity accuracy vs. parameter-count/datasize trade-off of
  /// broadcast-style ephemerides for lunar orbits.
  class CartesianEphemeris {
  public:
    CartesianEphemeris() = default;
    explicit CartesianEphemeris(EphemerisFitOptions options) : options_(options) {}

    /// @brief Fit this ephemeris model to a sampled trajectory.
    /// @param t_s  Sample epochs [s], strictly increasing, relative to any fixed
    ///             time origin.
    /// @param rv   Sampled Cartesian states [N x 6] ([m, m, m, m/s, m/s, m/s]),
    ///             row i corresponding to t_s(i).
    /// @return     Fitted parameter vector; see `ParamNames()` for the layout.
    VecXd Fit(const VecXd& t_s, const MatXd& rv) const;

    /// @brief Evaluate the fitted ephemeris at a set of epochs.
    /// @param t_s     Query epochs [s], same time origin as used in `Fit`.
    /// @param params  Parameter vector returned by `Fit`.
    /// @return        Evaluated Cartesian states [N x 6].
    MatXd Eval(const VecXd& t_s, const VecXd& params) const;

    /// @brief Evaluate the fit error of `params` against a reference trajectory.
    /// @param t_s     Reference epochs [s].
    /// @param rv_ref  Reference Cartesian states [N x 6] to compare against.
    /// @param params  Parameter vector returned by `Fit`.
    EphemerisFitErrorStats EvalError(const VecXd& t_s, const MatXd& rv_ref,
                                      const VecXd& params) const;

    /// @brief Number of scalar parameters in the vector returned by `Fit`.
    int NumParams() const;

    /// @brief Human-readable name of each parameter in the vector returned by
    /// `Fit` (same order/length), e.g. "t_ref", "t_fit", "a", "e", "i", "raan",
    /// "argp", "M_ref", "x_0", ..., "y_0", ..., "z_0", ...
    std::vector<std::string> ParamNames() const;

    const EphemerisFitOptions& GetOptions() const { return options_; }

  private:
    EphemerisFitOptions options_;
  };

}  // namespace lupnt
