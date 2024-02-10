#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <matplot/matplot.h>

namespace lupnt {

matplot::line_handle plot3(const VectorX &x, const VectorX &y, const VectorX &z,
                           std::string_view line_spec = "") {
  return matplot::plot3(EigenToStdVector(x), EigenToStdVector(y),
                        EigenToStdVector(z), line_spec);
}

matplot::line_handle plot(const VectorX &x, const VectorX &y,
                          std::string_view line_spec = "") {
  return matplot::plot(EigenToStdVector(x), EigenToStdVector(y), line_spec);
}

matplot::line_handle scatter3(const VectorX &x, const VectorX &y,
                              const VectorX &z, const VectorX &sizes,
                              std::string_view marker = "o") {
  return matplot::scatter3(EigenToStdVector(x), EigenToStdVector(y),
                           EigenToStdVector(z), EigenToStdVector(sizes),
                           marker);
}

}