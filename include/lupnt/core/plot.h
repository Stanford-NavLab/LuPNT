#include <matplot/matplot.h>

#include "lupnt/core/constants.h"

#pragma once

namespace lupnt {

  matplot::line_handle plot3(const VecX& x, const VecX& y, const VecX& z, double scale = 3,
                             std::string_view line_spec = "");
  matplot::line_handle plot3(const Vec3& xyz, std::string_view line_spec = "", double scale = 3);

  matplot::line_handle plot(const VecX& x, const VecX& y, std::string_view line_spec = "");

  matplot::line_handle scatter3(const VecX& x, const VecX& y, const VecX& z, const VecX& sizes,
                                double scale = 3, std::string_view marker = "o");

  matplot::surface_handle plot_body(NaifId body, Vec3 r_body = Vec3::Zero(), double scale = 3);

  void set_lim(Real lim, double scale = 3);
}  // namespace lupnt
