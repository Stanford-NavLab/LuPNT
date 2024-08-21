#include <matplot/matplot.h>

#include "lupnt/core/constants.h"

#pragma once

namespace lupnt {
  constexpr double default_scale = 3;

  matplot::line_handle Plot3(const VecX& x, const VecX& y, const VecX& z,
                             std::string_view line_spec = "", double scale = default_scale);
  matplot::line_handle Plot3(const Vec3& xyz, std::string_view line_spec = "",
                             double scale = default_scale);
  matplot::line_handle PlotArrow3(const Vec3& xyz, std::string_view line_spec = "",
                                  double scale = default_scale);
  matplot::line_handle Plot(const VecX& x, const VecX& y, std::string_view line_spec = "");

  matplot::surface_handle PlotBody(NaifId body, Vec3 r_body = Vec3::Zero(),
                                   double scale = default_scale);

  void SetLim(Real lim, double scale = default_scale);
}  // namespace lupnt
