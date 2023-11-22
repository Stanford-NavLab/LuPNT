#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) { return i + j; }
// void init_autodiff(py::module &m);
void init_coord_converter(py::module &m);
// // void init_constants(py::module &m);
void init_orbit_state(py::module &m);
// void init_dynamics(py::module &m);
// void init_math_utils(py::module &m);
// void init_spice_interface(py::module &m);

PYBIND11_MODULE(pylupnt, m) {
  // init_autodiff(m);
  init_coord_converter(m);
  // // init_constants(m);
  init_orbit_state(m);
  // init_math_utils(m);
  // init_spice_interface(m);
  // init_dynamics(m);
}
