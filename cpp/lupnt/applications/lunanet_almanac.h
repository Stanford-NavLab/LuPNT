#pragma once

#include <string>
#include <vector>

#include "lupnt/applications/lunanet_ephemeris.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Configuration for `Almanac::Fit`/`Eval`.
  struct AlmanacFitOptions {
    /// Degree of the secular polynomial fit to every element's time history
    /// (a, e, i, node, mean-anomaly residual, argument-of-latitude correction) in
    /// normalized time s = t_k / T_fit (beta_0 + beta_1 s + ... + beta_p s^p).
    /// Default 1 (linear) -- the argument-of-latitude parameterization (Algorithm 2
    /// of Iiyama & Gao) stays well-conditioned near periapsis even at linear order,
    /// so a linear fit is preferred to save message size. This is a coarse
    /// "GNSS-almanac-style" model, much lower accuracy than `CartesianEphemeris`,
    /// intended for a long validity window / small broadcast size rather than
    /// precision.
    int poly_order = 1;

    /// Number of Fourier harmonics added to each element's polynomial (a cos/sin
    /// pair per harmonic). The base angular frequency is element-specific -- the
    /// orbital mean motion (2*pi/T_orb) for the semi-major axis, and 4*pi/T_sid
    /// for the remaining (angular/eccentricity) elements -- capturing the periodic
    /// modulation that a low-order polynomial alone misses (per Iiyama & Gao,
    /// "Ephemeris and Almanac Design for Lunar Navigation Satellites"). Default 1.
    int num_fourier_terms = 1;

    /// Gravitational parameter of the central body [m^3/s^2].
    double gm = GM_MOON;

    /// Sidereal rotation period [s] of the central body, setting the base Fourier
    /// frequency 4*pi/sidereal_period_s for the non-semi-major-axis elements
    /// (default: the Moon's sidereal period, 27.321661 days).
    double sidereal_period_s = 27.321661 * SECS_DAY;

    /// Frame the input states (and the fitted output) are represented in. For an
    /// inertial frame (e.g. `Frame::MOON_CI`, the default) osculating elements are
    /// fit and reconstructed directly. For the rotating Moon-fixed principal-axis
    /// frame (`Frame::MOON_PA`) the osculating elements are computed in the
    /// Principal-Axis Inertial (PAI) frame (adding the omega x r velocity offset)
    /// and the reconstructed position/velocity is returned in MOON_PA (subtracting
    /// it back) -- per Iiyama & Gao, "Ephemeris and Almanac Design for Lunar
    /// Navigation Satellites". States passed to Fit/EvalError must already be
    /// expressed in this frame (convert e.g. from MOON_CI with `ConvertFrame`).
    Frame frame = Frame::MOON_CI;
  };

  /// @brief Coarse "almanac"-style orbit model (Algorithm 2 of Iiyama & Gao,
  /// "Ephemeris and Almanac Design for Lunar Navigation Satellites"): each
  /// osculating element (a, e, i, node) is a low-order polynomial plus an
  /// element-specific Fourier term (orbital mean motion for a, sidereal harmonic
  /// for the rest); the mean anomaly is the nominal two-body drift (integrated
  /// from the fitted semi-major axis) plus a polynomial/Fourier residual; and the
  /// argument of periapsis is replaced by an argument-of-latitude correction
  /// u = nu + (polynomial + Fourier), which keeps the reconstruction
  /// well-conditioned near periapsis of eccentric orbits. Position/velocity are
  /// built directly from (a, e, i, node, u) with velocity by finite differencing.
  /// This mirrors the coarse, long-validity almanac broadcast by GNSS
  /// constellations (e.g. GPS almanac pages), as distinct from the precise,
  /// short-validity `CartesianEphemeris`.
  ///
  /// Used directly (via the Python bindings) and by `EphemerisSimulation`
  /// (`lupnt/simulations/ephemeris/ephemeris_simulation.h`) to compare
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
    /// output): "t_ref", "t_fit", "a_ref", then per element the polynomial
    /// coefficients "a_p0", "a_p1", ... and Fourier coefficients "a_fc1", "a_fs1",
    /// ... for elements a, e, i, raan (node), M (mean-anomaly residual), and u
    /// (argument-of-latitude correction).
    std::vector<std::string> ParamNames() const;

    const AlmanacFitOptions& GetOptions() const { return options_; }

  private:
    AlmanacFitOptions options_;
  };

}  // namespace lupnt
