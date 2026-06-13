#include <matplot/matplot.h>

#include "lupnt/core/constants.h"

#pragma once

namespace lupnt {

  /// @brief Plot a 3D line through the points `(x, y, z)` using matplot++.
  ///
  /// Thin autodiff-to-double wrapper around `matplot::plot3`, used throughout
  /// `cpp/examples/` to visualize propagated trajectories (e.g. position
  /// columns of a `MatX6` orbit history) in 3D.
  ///
  /// @param x        X coordinates [consistent length units, e.g. km]
  /// @param y        Y coordinates [same units as `x`]
  /// @param z        Z coordinates [same units as `x`]
  /// @param line_spec matplot++ line/marker style string (e.g. `"b-"`, `"ro"`)
  /// @return         Handle to the created line, for further styling
  ///                 (`->display_name(...)`, `->line_width(...)`, etc.)
  matplot::line_handle Plot3(const VecX& x, const VecX& y, const VecX& z,
                             std::string_view line_spec = "");

  /// @brief `Plot3` overload taking an `[N x 3]` matrix of (x, y, z) columns.
  matplot::line_handle Plot3(const MatX& xyz, std::string_view line_spec);

  /// @brief Plot a single 3D point `xyz` as a marker (via `Plot3` with a
  /// single-element series).
  matplot::line_handle Scatter3(const Vec3& xyz, std::string_view line_spec = "");

  /// @brief Plot a single 2D point `xyz` as a marker (via `Plot`).
  matplot::line_handle Scatter(const Vec2& xyz, std::string_view line_spec = "");

  /// @brief `Scatter` overload taking an `[N x 2]` matrix of (x, y) columns,
  /// plotted as a 2D line/marker series.
  matplot::line_handle Scatter(const MatX2& xyz, std::string_view line_spec = "");

  /// @brief Plot a 3D arrow from `center` along direction `dir` (drawn as a
  /// 2-point line via `Plot3`).
  ///
  /// Used by `PlotFrame` to draw the three axes of a body/attitude frame for
  /// visual inspection of attitude/frame results in `cpp/examples/`.
  ///
  /// @param center   Arrow tail point [length units]
  /// @param dir      Arrow direction/length vector [same units as `center`]
  /// @param line_spec matplot++ line style string (e.g. `"r-"`)
  matplot::line_handle PlotArrow3(const Vec3& center, const Vec3& dir, std::string_view line_spec);

  /// @brief Draw the three axes of a rotation matrix `R` as colored arrows
  /// (red/green/blue for x/y/z) originating at `center`.
  ///
  /// Used in examples to visualize attitude/frame-rotation results, e.g.
  /// plotting a body or LVLH frame triad alongside a propagated trajectory.
  ///
  /// @param center Origin of the frame triad [length units]
  /// @param R      Rotation matrix whose rows are the frame's basis vectors
  ///               expressed in the plotting frame
  /// @return       Handles to the three drawn arrows (x, y, z axes)
  std::vector<matplot::line_handle> PlotFrame(const Vec3& center, const Mat3& R);

  /// @brief Plot a 2D line through the points `(x, y)` using matplot++.
  ///
  /// Thin autodiff-to-double wrapper around `matplot::plot`; used throughout
  /// `cpp/examples/` (e.g. time-history plots of states, residuals, errors).
  ///
  /// @param x        X coordinates (e.g. time)
  /// @param y        Y coordinates (e.g. state component, residual)
  /// @param line_spec matplot++ line/marker style string (e.g. `"b-"`, `"rx"`)
  /// @return         Handle to the created line, for further styling
  matplot::line_handle Plot(const VecX& x, const VecX& y, std::string_view line_spec = "");

  /// @brief Plot a celestial body as a shaded sphere centered at `r_body`.
  ///
  /// Used in scenario examples (e.g. `ex_frozen_orbits.cc`) to draw the Earth
  /// or Moon for visual context alongside plotted orbits, using the body's
  /// mean radius from `GetBodyRadius`.
  ///
  /// @param body   Body to draw (its radius is looked up via `GetBodyRadius`)
  /// @param r_body Center of the sphere [length units, e.g. 1e3 km to match
  ///               the axis labels set by this function]
  /// @return       Handle to the created surface plot
  matplot::surface_handle PlotBody(BodyId body, Vec3 r_body = Vec3::Zero());

  /// @brief Set symmetric `[-lim, lim]` limits on all three (x, y, z) axes of
  /// the current plot.
  void SetLim(Real lim);

  /// @brief Set `[lim_min, lim_max]` limits on all three (x, y, z) axes of the
  /// current plot.
  void SetLim(Real lim_min, Real lim_max);

}  // namespace lupnt
