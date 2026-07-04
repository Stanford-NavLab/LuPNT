#pragma once

#include <string>
#include <vector>

#include "lupnt/applications/ephemeris.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Configuration for `Almanac::Fit`/`Eval`.
  struct AlmanacFitOptions {
    /// Chebyshev polynomial order used for every orbital element's time history
    /// (a, e, i, raan, argp) and for the mean-anomaly residual left after removing
    /// the nominal two-body drift -- i.e. a coarse "GNSS-almanac-style" model,
    /// much lower order (and lower accuracy) than `CartesianEphemeris`, intended
    /// for a long validity window / small broadcast size rather than precision.
    int poly_order = 2;

    /// Gravitational parameter of the central body [m^3/s^2].
    double gm = GM_MOON;
  };

  /// @brief Coarse "almanac"-style orbit model: each classical orbital element
  /// (a, e, i, raan, argp) is represented by a low-order Chebyshev polynomial of
  /// time, and mean anomaly by the nominal two-body drift (integrated from the
  /// fitted, possibly time-varying, semi-major axis) plus a low-order polynomial
  /// residual. This mirrors the coarse, long-validity almanac broadcast by GNSS
  /// constellations (e.g. GPS almanac pages), as distinct from the precise,
  /// short-validity `CartesianEphemeris`.
  ///
  /// Used directly (via the Python bindings) and by `EphemerisSimulation`
  /// (`lupnt/simulations/Ephemeris/ephemeris_simulation.h`) to compare
  /// almanac-style vs. ephemeris-style broadcast data-size/accuracy trade-offs.
  class Almanac {
  public:
    Almanac() = default;
    explicit Almanac(AlmanacFitOptions options) : options_(options) {}

    /// @brief Fit this almanac model to a sampled trajectory.
    /// @param t_s  Sample epochs [s], strictly increasing.
    /// @param rv   Sampled Cartesian states [N x 6] ([m, m, m, m/s, m/s, m/s]).
    /// @return     Fitted parameter vector; see `ParamNames()` for the layout.
    VecXd Fit(const VecXd& t_s, const MatXd& rv) const;

    /// @brief Evaluate the fitted almanac at a set of epochs.
    MatXd Eval(const VecXd& t_s, const VecXd& params) const;

    /// @brief Evaluate the fit error of `params` against a reference trajectory.
    EphemerisFitErrorStats EvalError(const VecXd& t_s, const MatXd& rv_ref,
                                      const VecXd& params) const;

    /// @brief Number of scalar parameters in the vector returned by `Fit`.
    int NumParams() const;

    /// @brief Human-readable name of each parameter (same order as `Fit`'s
    /// output), e.g. "t_ref", "t_fit", "a_0", ..., "e_0", ..., "M_0", ...
    std::vector<std::string> ParamNames() const;

    const AlmanacFitOptions& GetOptions() const { return options_; }

  private:
    AlmanacFitOptions options_;
  };

}  // namespace lupnt
