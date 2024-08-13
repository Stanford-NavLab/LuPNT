#include "lupnt/core/plot.h"

#include <matplot/matplot.h>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {
  /// @brief Plot a 3D line
  /// @param x x-coordinates
  /// @param y y-coordinates
  /// @param z z-coordinates
  /// @param line_spec Line specification
  /// @return Line handle
  matplot::line_handle plot3(const VecX &x, const VecX &y, const VecX &z,
                             std::string_view line_spec) {
    return matplot::plot3(ToDouble(x), ToDouble(y), ToDouble(z), line_spec);
  }

  /// @brief Plot a 2D line
  /// @param x x-coordinates
  /// @param y y-coordinates
  /// @param line_spec Line specification
  /// @return Line handle
  matplot::line_handle plot(const VecX &x, const VecX &y, std::string_view line_spec) {
    return matplot::plot(ToDouble(x), ToDouble(y), line_spec);
  }

  /// @brief Plot a 3D scatter plot
  /// @param x x-coordinates
  /// @param y y-coordinates
  /// @param z z-coordinates
  /// @param sizes Marker sizes
  /// @param marker Marker style
  /// @return Line handle
  matplot::line_handle scatter3(const VecX &x, const VecX &y, const VecX &z, const VecX &sizes,
                                std::string_view marker) {
    return matplot::scatter3(ToDouble(x), ToDouble(y), ToDouble(z), ToDouble(sizes), marker);
  }
}  // namespace lupnt
