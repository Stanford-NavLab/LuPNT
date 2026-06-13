#include <matplot/matplot.h>

#include "lupnt/core/constants.h"
#include "lupnt/environment/body.h"
#include "lupnt/interfaces/matplot.h"

namespace lupnt {

  matplot::line_handle Plot3(const VecX &x, const VecX &y, const VecX &z,
                             std::string_view line_spec) {
    return matplot::plot3(x.cast<double>().eval(), y.cast<double>().eval(), z.cast<double>().eval(),
                          line_spec);
  }

  matplot::line_handle Plot3(const MatX &xyz, std::string_view line_spec) {
    VecX x = xyz.col(0), y = xyz.col(1), z = xyz.col(2);
    return Plot3(x, y, z, line_spec);
  }

  matplot::line_handle Scatter3(const Vec3 &xyz, std::string_view line_spec) {
    Vec1 x(xyz(0)), y(xyz(1)), z(xyz(2));
    return Plot3(x, y, z, line_spec);
  }

  matplot::line_handle Scatter(const Vec2 &xyz, std::string_view line_spec) {
    Vec1 x(xyz(0)), y(xyz(1));
    return Plot(x, y, line_spec);
  }
  matplot::line_handle Scatter(const MatX2 &xyz, std::string_view line_spec) {
    VecX x = xyz.col(0), y = xyz.col(1);
    return Plot(x, y, line_spec);
  }

  matplot::line_handle PlotArrow3(const Vec3 &center, const Vec3 &dir, std::string_view line_spec) {
    Vec2 x(center(0), center(0) + dir(0));
    Vec2 y(center(1), center(1) + dir(1));
    Vec2 z(center(2), center(2) + dir(2));
    return Plot3(x, y, z, line_spec);
  }

  matplot::line_handle Plot(const VecX &x, const VecX &y, std::string_view line_spec) {
    return matplot::plot(x.cast<double>().eval(), y.cast<double>().eval(), line_spec);
  }

  std::vector<matplot::line_handle> PlotFrame(const Vec3 &center, const Mat3 &R) {
    std::vector<std::string> line_specs = {"r-", "g-", "b-"};
    std::vector<matplot::line_handle> lines;
    for (int i = 0; i < 3; i++) lines.push_back(PlotArrow3(center, R.row(i), line_specs[i]));
    return lines;
  }

  matplot::surface_handle PlotBody(BodyId body, Vec3 r_body) {
    Vec3d r_body_ = r_body.cast<double>().eval();
    double radius = GetBodyRadius(body);
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
    h->edge_color("gray");
    xlabel("X [1e3 km]");
    ylabel("Y [1e3 km]");
    zlabel("Z [1e3 km]");
    // colormap(palette::gray());
    return h;
  }

  void SetLim(Real lim) {
    double lim_ = lim.val();
    matplot::xlim({-lim_, lim_});
    matplot::ylim({-lim_, lim_});
    matplot::zlim({-lim_, lim_});
  }

  void SetLim(Real lim_min, Real lim_max) {
    double lim_min_ = lim_min.val();
    double lim_max_ = lim_max.val();
    matplot::xlim({lim_min_, lim_max_});
    matplot::ylim({lim_min_, lim_max_});
    matplot::zlim({lim_min_, lim_max_});
  }

}  // namespace lupnt
