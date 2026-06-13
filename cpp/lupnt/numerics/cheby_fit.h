#pragma once
/**
 * @file cheby_fit.h
 * @brief Generic piecewise-Chebyshev polynomial model: fitting and evaluation.
 *
 * Extracted from the anonymous namespace of frame_conversions.cc so that
 * time_conversions.cc (TL−TT fit) and frame_conversions.cc (orientation fit)
 * both use the same infrastructure.
 *
 * ChebyshevFitSegment  — one polynomial segment in [t_mid±t_radius]
 * ChebyshevFitModel    — piecewise collection of segments, with Eval()
 * FitChebyshevModel()  — sample a function at Chebyshev-Gauss nodes and DCT-fit
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/interfaces/spice_cheby.h"  // cheby_eval_ad

namespace lupnt {

  /// @brief One Chebyshev-polynomial segment fitting an `num_dims`-dimensional function of
  /// time over [t_mid - t_radius, t_mid + t_radius].
  ///
  /// `coeffs[d]` holds the `num_coeffs` Chebyshev coefficients (T_n basis, constant term
  /// first) for dimension d -- directly usable by `cheby_eval_ad`. Produced by
  /// `FitChebyshevModel` and stored inside `ChebyshevFitModel::segments`.
  struct ChebyshevFitSegment {
    double t_mid = 0.0;
    double t_radius = 0.0;
    std::vector<std::vector<double>> coeffs;  // [dim][coeff index]
  };

  /// @brief A piecewise-Chebyshev fit of an `num_dims`-dimensional function of time over
  /// [t_start, t_end] (epochs in seconds relative to J2000).
  ///
  /// Built by `FitChebyshevModel` (or, for Euler angles, the analogous
  /// `FitEulerAngleModel` in frame_conversions.cc) and consulted via `Eval` whenever a
  /// SPICE-fitted quantity is requested -- e.g. by `InitFrameConversionFromSpice` for Earth
  /// orientation parameters and lunar libration angles, and by
  /// `GnssConstellation::GetSatelliteStateChebyshevEci` for fitted ECI ephemeris.
  struct ChebyshevFitModel {
    double t_start = 0.0;
    double t_end = 0.0;
    int num_dims = 0;
    int num_coeffs = 0;
    std::vector<ChebyshevFitSegment> segments;

    /// @brief Locate the segment of this model whose validity interval covers epoch `t`.
    ///
    /// Used internally by `Eval` to pick the correct set of Chebyshev coefficients before
    /// calling `cheby_eval_ad`.
    ///
    /// @param t  Epoch to look up [s, seconds relative to J2000]
    /// @return   Pointer to the covering segment, or `nullptr` if `t` lies outside
    ///           [t_start, t_end] or no segments are present (caller must not deref in
    ///           that case)
    const ChebyshevFitSegment* FindSegment(double t) const {
      if (segments.empty() || t < t_start || t > t_end) return nullptr;
      double span = t_end - t_start;
      int n = static_cast<int>(segments.size());
      int idx = span > 0.0 ? static_cast<int>((t - t_start) / span * n) : 0;
      idx = std::max(0, std::min(n - 1, idx));
      while (idx > 0 && t < segments[idx].t_mid - segments[idx].t_radius) idx--;
      while (idx + 1 < n && t > segments[idx].t_mid + segments[idx].t_radius) idx++;
      return &segments[idx];
    }

    /// @brief Evaluate the fitted `num_dims`-dimensional function (and optionally its exact
    /// analytic time derivative) at epoch `t`.
    ///
    /// Called by the frame-conversion routines (e.g. `RotPolarMotion`,
    /// `RotSideralMotion`/`RotSideralMotionDot`, `MoonCiToPa`/`MoonPaToCi`) once a
    /// SPICE-fitted model is active for the requested epoch, in place of the
    /// analytic/IERS-table-based computation.
    ///
    /// @param t   Epoch to evaluate at [s, seconds relative to J2000]
    /// @param f   Output: fitted function value, resized to `num_dims` [units of the
    ///            fitted quantity, e.g. rad for EOP/Euler angles, s for UT1-UTC]
    /// @param df  Optional output: exact analytic time derivative of `f` w.r.t. `t`, same
    ///            units per second; pass `nullptr` to skip
    /// @return    True if `t` is within [t_start, t_end] and `f` (and `df`, if requested)
    ///            were populated; false otherwise (outputs left unmodified)
    bool Eval(Real t, VecX* f, VecX* df) const {
      const ChebyshevFitSegment* seg = FindSegment(t.val());
      if (seg == nullptr) return false;
      double scale[2] = {seg->t_mid, seg->t_radius};
      f->resize(num_dims);
      if (df != nullptr) df->resize(num_dims);
      for (int d = 0; d < num_dims; d++) {
        Vec2 fd = cheby_eval_ad(t, scale, const_cast<double*>(seg->coeffs[d].data()), num_coeffs);
        (*f)(d) = fd(0);
        if (df != nullptr) (*df)(d) = fd(1);
      }
      return true;
    }
  };

  /// @brief Fit a piecewise Chebyshev polynomial model to a sampled `num_dims`-dimensional
  /// function of time.
  ///
  /// Samples `sample(t)` at the `num_coeffs` Chebyshev-Gauss nodes of [t_start, t_end],
  /// split into approximately `segment_length`-long equal segments, and fits a
  /// degree-(num_coeffs-1) Chebyshev polynomial to each segment via the standard discrete
  /// Chebyshev transform (DCT-II) -- i.e. the unique polynomial interpolant of `sample` at
  /// the Chebyshev-Gauss nodes, expressed in the T_n basis with the constant term first.
  ///
  /// Used by `InitFrameConversionFromSpice` to fit SPICE-derived Earth orientation
  /// parameters (via `ComputeEopFromSpice`) over a time window, and by
  /// `GnssConstellation`'s internal `FitStateHistoryChebyshev` helper to fit a satellite's
  /// ECI state history for fast ephemeris lookups (see
  /// `GnssConstellation::GetSatelliteStateChebyshevEci`).
  ///
  /// @param sample          Function returning the `num_dims`-dimensional sample vector at
  ///                         a given epoch [s, seconds relative to J2000]
  /// @param t_start         Start of the fit window [s, seconds relative to J2000]
  /// @param t_end           End of the fit window [s, seconds relative to J2000]
  /// @param num_dims        Dimensionality of `sample`'s output
  /// @param segment_length  Target length of each Chebyshev fit segment [s]
  /// @param num_coeffs      Number of Chebyshev coefficients (polynomial degree + 1) per
  ///                         dimension and segment
  /// @return                Fitted piecewise model, evaluable via `ChebyshevFitModel::Eval`
  inline ChebyshevFitModel FitChebyshevModel(const std::function<VecXd(double)>& sample,
                                             double t_start, double t_end, int num_dims,
                                             double segment_length, int num_coeffs) {
    if (!(t_end > t_start))
      throw std::runtime_error("FitChebyshevModel: t_end must be greater than t_start");
    if (num_coeffs < 2)
      throw std::runtime_error("FitChebyshevModel: num_coeffs must be at least 2");
    if (!(segment_length > 0.0))
      throw std::runtime_error("FitChebyshevModel: segment_length must be positive");

    double span = t_end - t_start;
    int num_segments = std::max(1, static_cast<int>(std::ceil(span / segment_length)));
    double seg_len = span / num_segments;

    ChebyshevFitModel model;
    model.t_start = t_start;
    model.t_end = t_end;
    model.num_dims = num_dims;
    model.num_coeffs = num_coeffs;
    model.segments.reserve(num_segments);

    for (int s = 0; s < num_segments; s++) {
      double seg_start = t_start + s * seg_len;
      double seg_end = (s == num_segments - 1) ? t_end : seg_start + seg_len;
      double mid = 0.5 * (seg_start + seg_end);
      double radius = 0.5 * (seg_end - seg_start);

      // Sample at the Chebyshev-Gauss nodes x_k = cos((2k+1)pi/(2N)),
      // mapped onto [seg_start, seg_end] via t = mid + radius*x_k.
      MatXd samples(num_coeffs, num_dims);
      for (int k = 0; k < num_coeffs; k++) {
        double theta = (2.0 * k + 1.0) * PI / (2.0 * num_coeffs);
        double t = mid + radius * std::cos(theta);
        samples.row(k) = sample(t).transpose();
      }

      // Discrete Chebyshev transform: c_0 = (1/N) sum_k f_k,
      // c_j = (2/N) sum_k f_k cos(j theta_k) for j > 0.
      ChebyshevFitSegment seg;
      seg.t_mid = mid;
      seg.t_radius = radius;
      seg.coeffs.assign(num_dims, std::vector<double>(num_coeffs, 0.0));
      for (int j = 0; j < num_coeffs; j++) {
        double w = (j == 0) ? 1.0 : 2.0;
        for (int k = 0; k < num_coeffs; k++) {
          double c = std::cos((2.0 * k + 1.0) * PI * j / (2.0 * num_coeffs));
          for (int d = 0; d < num_dims; d++) seg.coeffs[d][j] += samples(k, d) * c;
        }
        for (int d = 0; d < num_dims; d++) seg.coeffs[d][j] *= w / num_coeffs;
      }
      model.segments.push_back(std::move(seg));
    }
    return model;
  }

}  // namespace lupnt
