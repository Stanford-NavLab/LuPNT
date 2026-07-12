#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitFrameConverter(py::module &m) {
  py::enum_<Frame>(m, "Frame", "Reference frames (Earth, Moon, and solar-system planet frames).")
      .value("ITRF", Frame::ITRF, "International Terrestrial Reference Frame (Earth-fixed).")
      .value("ECEF", Frame::ECEF, "Earth-Centered Earth-Fixed (alias of ITRF).")
      .value("GCRF", Frame::GCRF, "Geocentric Celestial Reference Frame (Earth-centered inertial).")
      .value("ECI", Frame::ECI, "Earth-Centered Inertial (alias of EME/J2000).")
      .value("ICRF", Frame::ICRF,
             "International Celestial Reference Frame (solar-system-barycenter-centered).")
      .value("SER", Frame::SER, "Sun-Earth Rotating frame.")
      .value("GSE", Frame::GSE, "Geocentric Solar Ecliptic frame.")
      .value("EME", Frame::EME, "Earth-centered mean equator and equinox of J2000.")
      .value("MOD", Frame::MOD, "Mean-of-date equatorial frame.")
      .value("TOD", Frame::TOD, "True-of-date equatorial frame.")
      .value("EMR", Frame::EMR, "Earth-Moon Rotating frame.")
      .value("MOON_CI", Frame::MOON_CI, "Moon-centered inertial frame (ICRF-aligned axes).")
      .value("MOON_PA", Frame::MOON_PA, "Moon-fixed principal-axis frame.")
      .value("MOON_ME", Frame::MOON_ME, "Moon-fixed mean-Earth/polar-axis frame.")
      .value("MOON_OP", Frame::MOON_OP, "Moon orbit-plane frame.")
      // Solar-system planet body-fixed frames
      .value("MERCURY_FIXED", Frame::MERCURY_FIXED, "Mercury body-fixed frame.")
      .value("VENUS_FIXED", Frame::VENUS_FIXED, "Venus body-fixed frame.")
      .value("MARS_FIXED", Frame::MARS_FIXED, "Mars body-fixed frame.")
      .value("JUPITER_FIXED", Frame::JUPITER_FIXED, "Jupiter body-fixed frame.")
      .value("SATURN_FIXED", Frame::SATURN_FIXED, "Saturn body-fixed frame.")
      .value("URANUS_FIXED", Frame::URANUS_FIXED, "Uranus body-fixed frame.")
      .value("NEPTUNE_FIXED", Frame::NEPTUNE_FIXED, "Neptune body-fixed frame.")
      // Solar-system planet-centered inertial frames
      .value("MERCURY_CI", Frame::MERCURY_CI, "Mercury-centered inertial frame (ICRF-aligned).")
      .value("VENUS_CI", Frame::VENUS_CI, "Venus-centered inertial frame (ICRF-aligned).")
      .value("MARS_CI", Frame::MARS_CI, "Mars-centered inertial frame (ICRF-aligned).")
      .value("JUPITER_CI", Frame::JUPITER_CI, "Jupiter-centered inertial frame (ICRF-aligned).")
      .value("SATURN_CI", Frame::SATURN_CI, "Saturn-centered inertial frame (ICRF-aligned).")
      .value("URANUS_CI", Frame::URANUS_CI, "Uranus-centered inertial frame (ICRF-aligned).")
      .value("NEPTUNE_CI", Frame::NEPTUNE_CI, "Neptune-centered inertial frame (ICRF-aligned).")
      .export_values();

  // Vec6 = func(real, Vec6)
  m.def("convert_frame", py::overload_cast<Real, const Vec6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Convert a position+velocity state from `frame_in` to `frame_out` at `t_tdb` [m, m/s]. If "
        "`rotate_only`, skip the frame-origin translation.");
  // Vec3 = func(real, Vec3)
  m.def("convert_frame", py::overload_cast<Real, const Vec3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Convert a position from `frame_in` to `frame_out` at `t_tdb` [m]. If `rotate_only`, skip "
        "the frame-origin translation.");
  // VecX6 = func(VecX, Vec6)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const Vec6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Multi-epoch: convert the state `rv_in` from `frame_in` to `frame_out` at each epoch in "
        "`t_tdb` [N x 6].");
  // VecX3 = func(VecX, Vec3)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const Vec3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Multi-epoch: convert the position `r_in` from `frame_in` to `frame_out` at each epoch in "
        "`t_tdb` [N x 3].");
  // MatX6 = func(real, MatX6)
  m.def("convert_frame", py::overload_cast<Real, const MatX6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Convert each row (state) of `rv_in` from `frame_in` to `frame_out` at epoch `t_tdb` "
        "[N x 6].");
  // MatX3 = func(real, MatX3)
  m.def("convert_frame", py::overload_cast<Real, const MatX3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Convert each row (position) of `r_in` from `frame_in` to `frame_out` at epoch `t_tdb` "
        "[N x 3].");
  // MatX6 = func(VecX, MatX6)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const MatX6 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("rv_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Time-tagged rows: convert row i of `rv_in` from `frame_in` to `frame_out` at epoch "
        "t_tdb(i) [N x 6].");
  // MatX3 = func(VecX, MatX3)
  m.def("convert_frame",
        py::overload_cast<const VecX &, const MatX3 &, Frame, Frame, bool>(&ConvertFrame),
        py::arg("t_tdb"), py::arg("r_in"), py::arg("frame_in"), py::arg("frame_out"),
        py::arg("rotate_only") = false,
        "Time-tagged rows: convert row i of `r_in` from `frame_in` to `frame_out` at epoch "
        "t_tdb(i) [N x 3].");
  m.def("get_frame_rotation_translation",
        py::overload_cast<Real, Frame, Frame>(&GetFrameRotationTranslation), py::arg("t_tdb"),
        py::arg("frame_in"), py::arg("frame_out"),
        "Affine position transform (R, t) with r_out = R*r_in + t, from `frame_in` to `frame_out` "
        "at `t_tdb` (R dimensionless, t [m]).");
  m.def("get_frame_rotation_translation_rv",
        py::overload_cast<Real, Frame, Frame>(&GetFrameRotationTranslationRv), py::arg("t_tdb"),
        py::arg("frame_in"), py::arg("frame_out"),
        "Affine state transform (R6, t6) with x_out = R6*x_in + t6, from `frame_in` to "
        "`frame_out` at `t_tdb` (t6 [m, m/s]).");
  m.def("get_frame_center", &GetFrameCenter, py::arg("frame"),
        "NAIF body id of the frame's central body (e.g. GCRF/ITRF -> Earth, MOON_CI -> Moon, "
        "ICRF -> solar-system barycenter).");
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
