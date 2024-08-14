#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_frame_converter(py::module& m);
void init_constants(py::module& m);
void init_orbit_state(py::module& m);
void init_orbit_state_utils(py::module& m);
void init_math_utils(py::module& m);
void init_dynamics(py::module& m);
void init_measurements(py::module& m);
void init_spice_interface(py::module& m);

PYBIND11_MODULE(_pylupnt, m) {
  init_frame_converter(m);
  init_constants(m);
  init_orbit_state(m);
  init_orbit_state_utils(m);
  init_math_utils(m);
  init_dynamics(m);
  init_measurements(m);
  init_spice_interface(m);
}
