
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_coord_converter(py::module &m) {
  py::enum_<CoordSystem>(m, "CoordSystem")
      .value("ITRF", CoordSystem::ITRF)
      .value("GCRF", CoordSystem::GCRF)
      .value("ICRF", CoordSystem::ICRF)
      .value("SER", CoordSystem::SER)
      .value("GSE", CoordSystem::GSE)
      .value("EME", CoordSystem::EME)
      .value("MOD", CoordSystem::MOD)
      .value("TOD", CoordSystem::TOD)
      .value("EMR", CoordSystem::EMR)
      .value("MI", CoordSystem::MI)
      .value("PA", CoordSystem::PA)
      .value("ME", CoordSystem::ME)
      .value("OP", CoordSystem::OP)
      .export_values();

  py::class_<CoordConverter>(m, "CoordConverter")
      .def_static(
          "convert",
          [](double epoch, const Vector6d &rv_in, CoordSystem coord_sys_in,
             CoordSystem coord_sys_out) -> Vector6d {
            Vector6 rv_in_real = rv_in.cast<real>();
            return CoordConverter::Convert(epoch, rv_in_real, coord_sys_in,
                                           coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("rv_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      .def_static(
          "convert",
          [](double epoch, const Vector3d &r_in, CoordSystem coord_sys_in,
             CoordSystem coord_sys_out) -> Vector3d {
            Vector3 r_in_real = r_in.cast<real>();
            return CoordConverter::Convert(epoch, r_in_real, coord_sys_in,
                                           coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("r_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      .def_static(
          "convert",
          [](double epoch, const CartesianOrbitState &state_in,
             CoordSystem coord_sys_out) -> CartesianOrbitState {
            return CoordConverter::Convert(epoch, state_in, coord_sys_out);
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("state_in"),
          py::arg("coord_sys_out"));
}