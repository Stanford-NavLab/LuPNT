
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

MatXd def_convert_frame(double t_tai, const MatXd &rv_in, Frame frame_in, Frame frame_out) {
  Real t_tai_ = t_tai;
  if (rv_in.rows() == 3 && rv_in.cols() == 1) {  // 3x1
    Vec3 r_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, r_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() == 1 && rv_in.cols() == 3) {  // 1x3
    Vec3 r_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, r_in_, frame_in, frame_out).cast<double>().transpose();
  } else if (rv_in.rows() == 3 && rv_in.cols() > 1) {  // 3xN
    MatX3 rv_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() > 1 && rv_in.cols() == 3) {  // Nx3
    MatX3 rv_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>().transpose();
  } else if (rv_in.rows() == 6 && rv_in.cols() == 1) {  // 6x1
    Vec6 rv_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() == 1 && rv_in.cols() == 6) {  // 1x6
    Vec6 rv_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>().transpose();
  } else if (rv_in.rows() == 6 && rv_in.cols() > 1) {  // 6xN
    MatX6 rv_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() > 1 && rv_in.cols() == 6) {  // Nx6
    MatX6 rv_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>().transpose();
  }
  throw std::invalid_argument("Invalid input size for rv_in: " + std::to_string(rv_in.rows()) + "x"
                              + std::to_string(rv_in.cols()));
  return MatXd(0, 0);
}

MatXd def_convert_frame(VecXd t_tai, const MatXd &rv_in, Frame frame_in, Frame frame_out) {
  int n = t_tai.size();
  VecX t_tai_ = t_tai.cast<Real>().array();
  if (rv_in.rows() == 3 && rv_in.cols() == 1) {  // 3x1
    Vec3 r_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, r_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() == 1 && rv_in.cols() == 3) {  // 1x3
    Vec3 r_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, r_in_, frame_in, frame_out).cast<double>().transpose();
  } else if (rv_in.rows() == 3 && rv_in.cols() == n) {  // 3xN
    MatX3 rv_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() == n && rv_in.cols() == 3) {  // Nx3
    MatX3 rv_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>().transpose();
  } else if (rv_in.rows() == 6 && rv_in.cols() == 1) {  // 6x1
    Vec6 rv_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() == 1 && rv_in.cols() == 6) {  // 1x6
    Vec6 rv_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>().transpose();
  } else if (rv_in.rows() == 6 && rv_in.cols() == n) {  // 6xN
    MatX6 rv_in_ = rv_in.cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>();
  } else if (rv_in.rows() == n && rv_in.cols() == 6) {  // Nx6
    MatX6 rv_in_ = rv_in.transpose().cast<Real>();
    return ConvertFrame(t_tai_, rv_in_, frame_in, frame_out).cast<double>().transpose();
  }
  throw std::invalid_argument("Invalid input size for rv_in: " + std::to_string(rv_in.rows()) + "x"
                              + std::to_string(rv_in.cols()));
}

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

  // Vec6 = func(real, Vec6)
  m.def(
      "convert_frame",
      [](double t_tai, const Vec6d &rv_in, Frame frame_in, Frame frame_out) -> Vec6d {
        Vec6 rv_in_real = rv_in.cast<Real>();
        return ConvertFrame(t_tai, rv_in_real, frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("rv_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // Vec3 = func(real, Vec3)
  m.def(
      "convert_frame",
      [](double t_tai, const Vec3d &r_in, Frame frame_in, Frame frame_out) -> Vec3d {
        Vec3 r_in_real = r_in.cast<Real>();
        return ConvertFrame(t_tai, r_in_real, frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("r_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // CartesianOrbitState = func(real, CartesianOrbitState)
  m.def(
      "convert_frame",
      [](double t_tai, const CartesianOrbitState &state_in, Frame frame_out)
          -> CartesianOrbitState { return ConvertFrame(t_tai, state_in, frame_out); },
      "Convert frame", py::arg("t_tai"), py::arg("state_in"), py::arg("frame_out"));
  // Mat<-1,6> = func(VecX, Vec6)
  m.def(
      "convert_frame",
      [](const VecXd &t_tai, const Vec6d &rv_in, Frame frame_in, Frame frame_out) -> MatX6d {
        Vec6 rv_in_real = rv_in.cast<Real>();
        return ConvertFrame(t_tai, rv_in_real, frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("rv_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // Mat<-1,3> = func(VecX, Vec3)
  m.def(
      "convert_frame",
      [](const VecXd &t_tai, const Vec3d &r_in, Frame frame_in, Frame frame_out) -> MatX3d {
        Vec3 r_in_real = r_in.cast<Real>();
        return ConvertFrame(t_tai, r_in_real, frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("r_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // Mat<-1,6> = func(real, Mat<-1,6>)
  m.def(
      "convert_frame",
      [](double t_tai, const MatX6d &rv_in, Frame frame_in, Frame frame_out) -> MatX6d {
        return ConvertFrame(t_tai, rv_in.cast<Real>().eval(), frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("rv_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // Mat<-1,3> = func(real, Mat<-1,3>)
  m.def(
      "convert_frame",
      [](double t_tai, const MatX3d &r_in, Frame frame_in, Frame frame_out) -> MatX3d {
        return ConvertFrame(t_tai, r_in.cast<Real>().eval(), frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("r_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // Mat<-1,6> = func(VecX, Mat<-1,6>)
  m.def(
      "convert_frame",
      [](const VecXd &t_tai, const MatX6d &rv_in, Frame frame_in, Frame frame_out) -> MatX6d {
        return ConvertFrame(t_tai, rv_in.cast<Real>().eval(), frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("rv_in"), py::arg("frame_in"),
      py::arg("frame_out"));
  // Mat<-1,3> = func(VecX, Mat<-1,3>)
  m.def(
      "convert_frame",
      [](const VecXd &t_tai, const MatX3d &r_in, Frame frame_in, Frame frame_out) -> MatX3d {
        return ConvertFrame(t_tai, r_in.cast<Real>().eval(), frame_in, frame_out).cast<double>();
      },
      "Convert frame", py::arg("t_tai"), py::arg("r_in"), py::arg("frame_in"),
      py::arg("frame_out"));
}
