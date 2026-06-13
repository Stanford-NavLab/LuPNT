#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <lupnt/lupnt.h>

#include <any>
#include <functional>
#include <iostream>
#include <vector>

#include "py_pybind11.h"
// #include "py_real.h"
#include "py_eigen.h"

using namespace lupnt;

py::tuple JacobianWrapper(py::function func_py, const py::tuple& wrts_py, const py::tuple& ats_py) {
  VecX y;
  std::vector<MatX> Js;
  Js.reserve(wrts_py.size());

  // Iterate over wrt
  for (ssize_t i = 0; i < wrts_py.size(); ++i) {
    VecX x = wrts_py[i].cast<VecX>();

    // Find in at
    bool found = false;
    ssize_t j = 0;
    for (; j < ats_py.size(); ++j) {
      if (ats_py[j].cast<VecX>() == x) {
        found = true;
        break;
      }
    }
    if (!found) throw py::value_error("wrt variable not found in ats");

    // Callable
    auto func = [&](const VecX& xx) -> VecX {
      py::tuple args(ats_py.size());
      args[j] = py::cast(xx, py::return_value_policy::reference);
      for (ssize_t k = 0; k < ats_py.size(); ++k) {
        if (k != j) args[k] = ats_py[k];
      }
      py::object out = func_py(*args);

      if (py::isinstance<py::float_>(out) || py::isinstance<py::int_>(out)
          || py::isinstance<Real>(out)) {
        VecX v(1);
        v << out.cast<Real>();
        return v;
      }
      return out.cast<VecX>();
    };

    MatX J;
    autodiff::jacobian(func, autodiff::wrt(x), autodiff::at(x), y, J);
    Js.push_back(J);
  }

  // Return
  if (y.size() == 1) {
    return py::make_tuple(y(0), Js);
  }
  return py::make_tuple(y, Js);
}

void InitAutodiff(py::module& m) {
  // m.def("jacobian", &JacobianWrapper, py::arg("func"), py::arg("wrt"), py::arg("at"));
  // ExportReal<1, double>(m, "Real");
  // Add implicit conversions from Python numeric types to Real
  // py::implicitly_convertible<py::float_, Real>();
  // py::implicitly_convertible<py::int_, Real>();
}
