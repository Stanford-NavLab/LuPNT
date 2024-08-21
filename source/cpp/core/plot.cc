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
  matplot::line_handle plot3(const VecX &x, const VecX &y, const VecX &z, double scale,
                             std::string_view line_spec) {
    scale = pow(10, scale);
    return matplot::plot3(ToDouble(x / scale), ToDouble(y / scale), ToDouble(z / scale), line_spec);
  }

  matplot::line_handle plot3(const Vec3 &xyz, std::string_view line_spec, double scale) {
    scale = pow(10, scale);
    Vec2d x(0, xyz(0).val() / scale), y(0, xyz(1).val() / scale), z(0, xyz(2).val() / scale);
    return matplot::plot3(x, y, z, line_spec);
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
                                double scale, std::string_view marker) {
    return matplot::scatter3(ToDouble(x / scale), ToDouble(y / scale), ToDouble(z / scale),
                             ToDouble(sizes), marker);
  }

  matplot::surface_handle plot_body(NaifId body, Vec3 r_body, double scale) {
    scale = pow(10, scale);
    Vec3d r_body_ = ToDouble(r_body / scale);
    double radius = GetBodyRadius(body) / scale;
    using namespace matplot;
    int n = 40;
    auto theta = linspace(0, pi, n);
    auto phi = linspace(0, 2 * pi, n);
    auto [T, P] = meshgrid(theta, phi);
    auto X = transform(T, P, [radius, r_body_](double theta, double phi) {
      return radius * sin(theta) * cos(phi) + r_body_(0);
    });
    auto Y = transform(T, P, [radius, r_body_](double theta, double phi) {
      return radius * sin(theta) * sin(phi) + r_body_(1);
    });
    auto Z = transform(T, P, [radius, r_body_](double theta, double phi) {
      (void)phi;
      return radius * cos(theta) + r_body_(2);
    });
    auto h = surf(X, Y, Z);
    h->edge_color("none");
    xlabel("X [1e3 km]");
    ylabel("Y [1e3 km]");
    zlabel("Z [1e3 km]");
    // colormap(palette::gray());
    return h;
  }

  void set_lim(Real lim, double scale) {
    scale = pow(10, scale);
    double lim_ = lim.val() / scale;
    matplot::xlim({-lim_, lim_});
    matplot::ylim({-lim_, lim_});
    matplot::zlim({-lim_, lim_});
  }

}  // namespace lupnt
