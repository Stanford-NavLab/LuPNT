
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_coord_converter(py::module &m) {
  py::enum_<CoordSystem>(m, "CoordSystem")
      .value("ITRF", CoordSystem::ITRF)
      .value("ECEF", CoordSystem::ECEF)
      .value("GCRF", CoordSystem::GCRF)
      .value("ECI", CoordSystem::ECI)
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
      // Vector6 = func(real, Vector6)
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
      // Vector3 = func(real, Vector3)
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
      // CartesianOrbitState = func(real, CartesianOrbitState)
      .def_static(
          "convert",
          [](double epoch, const CartesianOrbitState &state_in,
             CoordSystem coord_sys_out) -> CartesianOrbitState {
            return CoordConverter::Convert(epoch, state_in, coord_sys_out);
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("state_in"),
          py::arg("coord_sys_out"))
      // Matrix<-1,6> = func(VectorX, Vector6)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Vector6d &rv_in,
             CoordSystem coord_sys_in, CoordSystem coord_sys_out) -> MatrixXd {
            Vector6 rv_in_real = rv_in.cast<real>();
            return CoordConverter::Convert(epoch, rv_in_real, coord_sys_in,
                                           coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("rv_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      // Matrix<-1,3> = func(VectorX, Vector3)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Vector3d &r_in,
             CoordSystem coord_sys_in, CoordSystem coord_sys_out) -> MatrixXd {
            Vector3 r_in_real = r_in.cast<real>();
            return CoordConverter::Convert(epoch, r_in_real, coord_sys_in,
                                           coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("r_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      // Matrix<-1,6> = func(real, Matrix<-1,6>)
      .def_static(
          "convert",
          [](double epoch, const Matrixd<-1, 6> &rv_in,
             CoordSystem coord_sys_in,
             CoordSystem coord_sys_out) -> Matrixd<-1, 6> {
            return CoordConverter::Convert(epoch, rv_in.cast<real>().eval(),
                                           coord_sys_in, coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("rv_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      // Matrix<-1,3> = func(real, Matrix<-1,3>)
      .def_static(
          "convert",
          [](double epoch, const Matrixd<-1, 3> &r_in, CoordSystem coord_sys_in,
             CoordSystem coord_sys_out) -> Matrixd<-1, 3> {
            return CoordConverter::Convert(epoch, r_in.cast<real>().eval(),
                                           coord_sys_in, coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("r_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      // Matrix<-1,6> = func(VectorX, Matrix<-1,6>)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Matrixd<-1, 6> &rv_in,
             CoordSystem coord_sys_in,
             CoordSystem coord_sys_out) -> Matrixd<-1, 6> {
            return CoordConverter::Convert(epoch, rv_in.cast<real>().eval(),
                                           coord_sys_in, coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("rv_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"))
      // Matrix<-1,3> = func(VectorX, Matrix<-1,3>)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Matrixd<-1, 3> &r_in,
             CoordSystem coord_sys_in,
             CoordSystem coord_sys_out) -> Matrixd<-1, 3> {
            return CoordConverter::Convert(epoch, r_in.cast<real>().eval(),
                                           coord_sys_in, coord_sys_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("epoch"), py::arg("r_in"),
          py::arg("coord_sys_in"), py::arg("coord_sys_out"));
}