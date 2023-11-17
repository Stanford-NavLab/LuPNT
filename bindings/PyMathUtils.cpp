#include <lupnt/numerics/MathUtils.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace LPT;

void init_math_utils(py::module &m) {
  // m.def("wrapToPi", &wrapToPi, "Wrap angle to [-pi, pi]");
  m.def(
      "wrapToPi", [](double angle) -> double { return wrapToPi(angle); },
      "Wrap angle to [-pi, pi]");
}
