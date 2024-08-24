#include <pybind11/pybind11.h>

#include "py_constants.cc"
#include "py_dynamics.cc"
#include "py_forces.cc"
#include "py_frame_converter.cc"
#include "py_kernels.cc"
#include "py_math_utils.cc"
#include "py_measurements.cc"
#include "py_orbit_state.cc"
#include "py_orbit_state_utils.cc"
#include "py_sandbox.cc"
#include "py_spice_interface.cc"
#include "py_time_converter.cc"

namespace py = pybind11;

PYBIND11_MODULE(_pylupnt, m) {
  init_frame_converter(m);
  init_constants(m);
  init_orbit_state(m);
  init_orbit_state_utils(m);
  init_math_utils(m);
  init_dynamics(m);
  init_measurements(m);
  init_spice_interface(m);
  init_time_converter(m);
  init_kernels(m);
  init_sandbox(m);
  init_forces(m);
}
