#include <lupnt/numerics/math_utils.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_math_utils(py::module &m) {
  // m.def("wrapToPi", &wrapToPi, "Wrap angle to [-pi, pi]");
  m.def(
      "wrapToPi", [](double angle) -> double { return wrapToPi(angle).val(); },
      "Wrap angle to [-pi, pi]");
  m.def(
      "wrapTo2Pi",
      [](double angle) -> double { return wrapTo2Pi(angle).val(); },
      "Wrap angle to [0, 2pi]");
  m.def(
      "wrapToPi",
      [](const VectorXd &angle) -> VectorXd {
        return wrapToPi(angle).cast<double>();
      },
      "Wrap angle to [-pi, pi]");
  m.def(
      "wrapTo2Pi",
      [](const VectorXd &angle) -> VectorXd {
        return wrapTo2Pi(angle).cast<double>();
      },
      "Wrap angle to [0, 2pi]");
}
