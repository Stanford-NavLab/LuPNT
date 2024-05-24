// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/physics/frame_converter.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_frame_converter(py::module &m) {
  py::enum_<Frame>(m, "Frame")
      .value("ITRF", Frame::ITRF)
      .value("GCRF", Frame::GCRF)
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
      .export_values();

  py::class_<FrameConverter>(m, "FrameConverter")
      .def_static(
          "convert",
          [](double epoch, const Vector6d &rv_in, Frame frame_in,
             Frame frame_out) -> Vector6d {
            return FrameConverter::Convert(epoch, rv_in, frame_in, frame_out)
                .cast<double>();
          },
          "Convert coordinate system", py::arg("rv_in"), py::arg("epoch"),
          py::arg("frame_in"), py::arg("frame_out"));
}