
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
      .value("MOON_CI", Frame::MOON_CI)
      .value("MOON_PA", Frame::MOON_PA)
      .value("MOON_ME", Frame::MOON_ME)
      .value("MOON_OP", Frame::MOON_OP)
      .value("MARS_FIXED", Frame::MARS_FIXED)
      .value("VENUS_FIXED", Frame::VENUS_FIXED)
      .export_values();

  py::class_<FrameConverter>(m, "FrameConverter")
      // Vec6 = func(real, Vec6)
      .def_static(
          "convert",
          [](double epoch, const Vec6d &rv_in, Frame frame_in,
             Frame frame_out) -> Vec6d {
            Vec6 rv_in_real = rv_in.cast<real>();
            return FrameConverter::Convert(epoch, rv_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Vec3 = func(real, Vec3)
      .def_static(
          "convert",
          [](double epoch, const Vec3d &r_in, Frame frame_in,
             Frame frame_out) -> Vec3d {
            Vec3 r_in_real = r_in.cast<real>();
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
      // Mat<-1,6> = func(VecX, Vec6)
      .def_static(
          "convert",
          [](const VecXd &epoch, const Vec6d &rv_in, Frame frame_in,
             Frame frame_out) -> VecXd {
            Vec6 rv_in_real = rv_in.cast<real>();
            return FrameConverter::Convert(epoch, rv_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Mat<-1,3> = func(VecX, Vec3)
      .def_static(
          "convert",
          [](const VecXd &epoch, const Vec3d &r_in, Frame frame_in,
             Frame frame_out) -> VecXd {
            Vec3 r_in_real = r_in.cast<real>();
            return FrameConverter::Convert(epoch, r_in_real, frame_in,
                                           frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Mat<-1,6> = func(real, Mat<-1,6>)
      .def_static(
          "convert",
          [](double epoch, const Matd<-1, 6> &rv_in, Frame frame_in,
             Frame frame_out) -> Matd<-1, 6> {
            return FrameConverter::Convert(epoch, rv_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Mat<-1,3> = func(real, Mat<-1,3>)
      .def_static(
          "convert",
          [](double epoch, const Matd<-1, 3> &r_in, Frame frame_in,
             Frame frame_out) -> Matd<-1, 3> {
            return FrameConverter::Convert(epoch, r_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Mat<-1,6> = func(VecX, Mat<-1,6>)
      .def_static(
          "convert",
          [](const VecXd &epoch, const Matd<-1, 6> &rv_in, Frame frame_in,
             Frame frame_out) -> Matd<-1, 6> {
            return FrameConverter::Convert(epoch, rv_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("rv_in"),
          py::arg("frame_in"), py::arg("frame_out"))
      // Mat<-1,3> = func(VecX, Mat<-1,3>)
      .def_static(
          "convert",
          [](const VecXd &epoch, const Matd<-1, 3> &r_in, Frame frame_in,
             Frame frame_out) -> Matd<-1, 3> {
            return FrameConverter::Convert(epoch, r_in.cast<real>().eval(),
                                           frame_in, frame_out)
                .cast<double>();
          },
          "Convert frame", py::arg("epoch"), py::arg("r_in"),
          py::arg("frame_in"), py::arg("frame_out"));
}