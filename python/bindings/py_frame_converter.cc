#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitFrameConverter(py::module &m) {
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
      // Solar-system planet body-fixed frames
      .value("MERCURY_FIXED", Frame::MERCURY_FIXED)
      .value("VENUS_FIXED", Frame::VENUS_FIXED)
      .value("MARS_FIXED", Frame::MARS_FIXED)
      .value("JUPITER_FIXED", Frame::JUPITER_FIXED)
      .value("SATURN_FIXED", Frame::SATURN_FIXED)
      .value("URANUS_FIXED", Frame::URANUS_FIXED)
      .value("NEPTUNE_FIXED", Frame::NEPTUNE_FIXED)
      // Solar-system planet-centered inertial frames
      .value("MERCURY_CI", Frame::MERCURY_CI)
      .value("VENUS_CI", Frame::VENUS_CI)
      .value("MARS_CI", Frame::MARS_CI)
      .value("JUPITER_CI", Frame::JUPITER_CI)
      .value("SATURN_CI", Frame::SATURN_CI)
      .value("URANUS_CI", Frame::URANUS_CI)
      .value("NEPTUNE_CI", Frame::NEPTUNE_CI)
      .export_values();

  // Vec6 = func(real, Vec6)
  m.def("convert_frame", py::overload_cast<Real, const Vec6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // Vec3 = func(real, Vec3)
  m.def("convert_frame", py::overload_cast<Real, const Vec3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // VecX6 = func(VecX, Vec6)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const Vec6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // VecX3 = func(VecX, Vec3)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const Vec3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // MatX6 = func(real, MatX6)
  m.def("convert_frame", py::overload_cast<Real, const MatX6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // MatX3 = func(real, MatX3)
  m.def("convert_frame", py::overload_cast<Real, const MatX3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // MatX6 = func(VecX, MatX6)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const MatX6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  // MatX3 = func(VecX, MatX3)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const MatX3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false);
  m.def("get_frame_rotation_translation",
        py::overload_cast<Real, Frame, Frame>(&GetFrameRotationTranslation), py::arg("t_tdb"),
        py::arg("frame_in"), py::arg("frame_out"));
  m.def("get_frame_rotation_translation_rv",
        py::overload_cast<Real, Frame, Frame>(&GetFrameRotationTranslationRv), py::arg("t_tdb"),
        py::arg("frame_in"), py::arg("frame_out"));
  m.def("get_frame_center", &GetFrameCenter, py::arg("frame"));
  m.def(
      "compute_eop_from_spice",
      [](double t_tdb) -> Vec3d {
        SpiceEopParams p = ComputeEopFromSpice(Real(t_tdb));
        return Vec3d(p.x_pole.val(), p.y_pole.val(), p.ut1_utc.val());
      },
      "Derive (x_pole [rad], y_pole [rad], UT1-UTC [s]) at t_tdb from SPICE", py::arg("t_tdb"));
  m.def(
      "get_lunar_orientation_angles",
      [](double t_tdb) -> Vec6d {
        Real t_tdb_ = Real(t_tdb);
        Vec6 angles = GetLunarMantleData(t_tdb_);
        return angles.cast<double>();
      },
      "Get lunar mantle rotation angles", py::arg("t_tdb"));
  //   m.def(
  //       "get_lunar_orientation_angles",
  //       [](double t_tdb) -> Vec6d {
  //         Real t_tdb_ = Real(t_tdb);
  //         Vec6 angles = GetLunarMantleData(t_tdb_);
  //         return angles.cast<double>();
  //       },
  //       "Get lunar mantle rotation angles", py::arg("t_tdb"));
}
