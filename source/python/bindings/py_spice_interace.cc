#include <lupnt/core/constants.h>
#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/spice_interface.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module &m) {
  py::class_<SpiceInterface>(m, "SpiceInterface")
      .def_static("load_spice_kernel", &lupnt::spice::LoadSpiceKernel)
      .def_static("extract_pck_coeffs", &lupnt::spice::ExtractPckCoeffs)
      .def_static(
          "get_frame_conversion_mat",
          [](double t_tai, lupnt::Frame from, lupnt::Frame to) -> lupnt::Mat6d {
            return lupnt::spice::GetFrameConversionMat(t_tai, from, to).cast<double>();
          },
          py::arg("t_tai"), py::arg("from"), py::arg("to"))
      .def_static(
          "string_to_tdb",
          [](std::string utc) -> double { return lupnt::spice::String2TDB(utc).val(); },
          py::arg("gregorian_date"))
      .def_static(
          "string_to_tai",
          [](std::string utc) -> double { return lupnt::spice::String2TAI(utc).val(); },
          py::arg("gregorian_date"))
      .def_static(
          "tdb_to_string_utc",
          [](double t_tdb, int prec) { return lupnt::spice::TDBtoStringUTC(t_tdb, prec); },
          py::arg("t_tdb"), py::arg("precision"))
      .def_static(
          "convert_time",
          [](double t_tai, std::string from, std::string to) {
            return lupnt::spice::ConvertTime(t_tai, from, to).val();
          },
          py::arg("t_tai"), py::arg("from"), py::arg("to"))
      .def_static(
          "get_body_pos_vel",
          [](double t_tai, lupnt::NaifId center, lupnt::NaifId target) -> lupnt::Vec6d {
            return lupnt::spice::GetBodyPosVel(t_tai, center, target).cast<double>();
          },
          py::arg("t_tai"), py::arg("center"), py::arg("target"))
      .def_static(
          "get_body_pos_vel",
          [](lupnt::VecXd t_tai, lupnt::NaifId center, lupnt::NaifId target) -> lupnt::Matd<-1, 6> {
            return lupnt::spice::GetBodyPosVel(t_tai, center, target).cast<double>();
          },
          py::arg("t_tai"), py::arg("center"), py::arg("target"))
      .def_static(
          "get_body_pos_spice",
          [](double t_tai, lupnt::NaifId obs, lupnt::NaifId target, lupnt::Frame ref_frame,
             std::string ab_orrection) -> lupnt::Vec3d {
            return lupnt::spice::GetBodyPosSpice(t_tai, obs, target, ref_frame, ab_orrection)
                .cast<double>();
          },
          py::arg("t_tai"), py::arg("obs"), py::arg("target"), py::arg("ref_frame"),
          py::arg("ab_correction"));
}