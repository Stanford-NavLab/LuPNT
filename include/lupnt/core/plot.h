#include <matplot/matplot.h>

#include "lupnt/core/constants.h"

#pragma once

namespace lupnt {

  matplot::line_handle plot3(const VecX& x, const VecX& y, const VecX& z,
                             std::string_view line_spec = "");

  matplot::line_handle plot(const VecX& x, const VecX& y, std::string_view line_spec = "");

  matplot::line_handle scatter3(const VecX& x, const VecX& y, const VecX& z, const VecX& sizes,
                                std::string_view marker = "o");
}  // namespace lupnt
