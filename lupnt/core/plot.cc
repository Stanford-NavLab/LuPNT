#include "lupnt/core/plot.h"

#include <matplot/matplot.h>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {
matplot::line_handle plot3(const VecX &x, const VecX &y, const VecX &z,
                           std::string_view line_spec) {
  return matplot::plot3(EigenToStdVec(x), EigenToStdVec(y), EigenToStdVec(z),
                        line_spec);
}

matplot::line_handle plot(const VecX &x, const VecX &y,
                          std::string_view line_spec) {
  return matplot::plot(EigenToStdVec(x), EigenToStdVec(y), line_spec);
}

matplot::line_handle scatter3(const VecX &x, const VecX &y, const VecX &z,
                              const VecX &sizes, std::string_view marker) {
  return matplot::scatter3(EigenToStdVec(x), EigenToStdVec(y), EigenToStdVec(z),
                           EigenToStdVec(sizes), marker);
}
}  // namespace lupnt