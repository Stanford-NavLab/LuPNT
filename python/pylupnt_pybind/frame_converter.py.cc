
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_frame_converter(py::module &m) {
  py::enum_<Frame>(m, "Frame")
      .value("ITRF", Frame::ITRF)
      .value("ECEF", Frame::ECEF)
      .value("GCRF", Frame::GCRF)
      .value("ECI", Frame::ECI)
      .value("ICRF", Frame::ICRF)
      .value("SER", Frame::SER)
      .value("GSE", Frame::GSE)
      .value("EME", Frame::EME)
      .value("MOD", Frame::MOD)
      .value("TOD", Frame::TOD)
      .value("EMR", Frame::EMR)
      .value("MI", Frame::MI)
      .value("PA", Frame::PA)
      .value("ME", Frame::ME)
      .value("OP", Frame::OP)
      .export_values();

  py::class_<FrameConverter>(m, "FrameConverter")
      // Vector6 = func(real, Vector6)
      .def_static(
          "convert",
          [](double epoch, const Vector6d &rv_in, Frame frame_in,
             Frame frame_out) -> Vector6d {
            Vector6 rv_in_real = rv_in.cast<real>();
            return FrameConverter::Convert(epoch, rv_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Vector3 = func(real, Vector3)
      .def_static(
          "convert",
          [](double epoch, const Vector3d &r_in, Frame frame_in,
             Frame frame_out) -> Vector3d {
            Vector3 r_in_real = r_in.cast<real>();
            return FrameConverter::Convert(epoch, r_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // CartesianOrbitState = func(real, CartesianOrbitState)
      .def_static(
          "convert",
          [](double epoch, const CartesianOrbitState &state_in,
             Frame frame_out) -> CartesianOrbitState {
            return FrameConverter::Convert(epoch, state_in, frame_out);
          },
          "Convert frame", py::arg("epoch"), py::arg("state_in"),
          py::arg("frame_out"))
      // Matrix<-1,6> = func(VectorX, Vector6)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Vector6d &rv_in, Frame frame_in,
             Frame frame_out) -> MatrixXd {
            Vector6 rv_in_real = rv_in.cast<real>();
            return FrameConverter::Convert(epoch, rv_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Matrix<-1,3> = func(VectorX, Vector3)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Vector3d &r_in, Frame frame_in,
             Frame frame_out) -> MatrixXd {
            Vector3 r_in_real = r_in.cast<real>();
            return FrameConverter::Convert(epoch, r_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Matrix<-1,6> = func(real, Matrix<-1,6>)
      .def_static(
          "convert",
          [](double epoch, const Matrixd<-1, 6> &rv_in, Frame frame_in,
             Frame frame_out) -> Matrixd<-1, 6> {
            return FrameConverter::Convert(epoch, rv_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Matrix<-1,3> = func(real, Matrix<-1,3>)
      .def_static(
          "convert",
          [](double epoch, const Matrixd<-1, 3> &r_in, Frame frame_in,
             Frame frame_out) -> Matrixd<-1, 3> {
            return FrameConverter::Convert(epoch, r_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Matrix<-1,6> = func(VectorX, Matrix<-1,6>)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Matrixd<-1, 6> &rv_in, Frame frame_in,
             Frame frame_out) -> Matrixd<-1, 6> {
            return FrameConverter::Convert(epoch, rv_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Matrix<-1,3> = func(VectorX, Matrix<-1,3>)
      .def_static(
          "convert",
          [](const VectorXd &epoch, const Matrixd<-1, 3> &r_in, Frame frame_in,
             Frame frame_out) -> Matrixd<-1, 3> {
            return FrameConverter::Convert(epoch, r_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"));
}