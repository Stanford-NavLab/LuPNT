#include <lupnt/core/constants.h>
#include <matplot/matplot.h>

#pragma once

namespace lupnt {
matplot::line_handle plot3(const VectorX &x, const VectorX &y, const VectorX &z,
                           std::string_view line_spec = "");

matplot::line_handle plot(const VectorX &x, const VectorX &y,
                          std::string_view line_spec = "");

matplot::line_handle scatter3(const VectorX &x, const VectorX &y,
                              const VectorX &z, const VectorX &sizes,
                              std::string_view marker = "o");
}  // namespace lupnt
