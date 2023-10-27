#include "pybind11.hxx"

// dynamics
void init_autodiff(py::module &m);
void init_coord_converter(py::module &m);
void init_spice_interface(py::module &m);
void init_state(py::module &m);
void init_constants(py::module &m);
void init_dynamics(py::module &m);
void init_math_utils(py::module &m);

PYBIND11_MODULE(pylupnt, m) {
  // dynamics
  init_autodiff(m);         // PyAutodiff.cpp
  init_math_utils(m);       // PyMathUtils.cpp
  init_constants(m);        // PyConstants.cpp
  init_coord_converter(m);  // PyCoordConverter.cpp
  init_spice_interface(m);  // PySpiceInterface.cpp
  init_state(m);            // PyState.cpp
  init_dynamics(m);         // PyDynamics.cpp
}