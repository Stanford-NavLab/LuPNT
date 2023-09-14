#include <lupnt/numerics/MathUtils.h>
#include <pybind11/pybind11.h>

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace py = pybind11;
using namespace LPT;

void init_math_utils(py::module &m) {
  m.def("wrapToPi", &wrapToPi, "Wrap angle to [-pi, pi]");
}
