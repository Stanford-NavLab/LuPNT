#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module& m) {
  auto m_spice = m.def_submodule("spice", "This is A.");

  m_spice
      .def("load_spice_kernel",
           static_cast<void (*)(const std::filesystem::path&)>(&lupnt::spice::LoadSpiceKernel))
      .def("extract_pck_coeffs", &lupnt::spice::ExtractPckCoeffs)
      .def(
          "get_frame_conversion_mat",
          [](double t_tdb, const std::string& from, const std::string& to) -> lupnt::Mat6d {
            return lupnt::spice::GetFrameConversionMat(t_tdb, from, to).cast<double>();
          },
          py::arg("t_tdb"), py::arg("from"), py::arg("to"))
      .def(
          "string2tdb",
          [](std::string utc) -> double { return lupnt::spice::StringToTdb(utc).val(); },
          py::arg("gregorian_date"))
      .def(
          "string2tai",
          [](std::string utc) -> double { return lupnt::spice::StringToTai(utc).val(); },
          py::arg("gregorian_date"))
      .def(
          "tdb2string_utc",
          [](double t_tdb, int prec) { return lupnt::spice::TDBtoStringUTC(t_tdb, prec); },
          py::arg("t_tdb"), py::arg("precision"))
      .def("convert_time", &lupnt::spice::ConvertTime, py::arg("t"), py::arg("from"), py::arg("to"))
      .def(
          "get_body_pos_vel",
          [](double t_tdb, lupnt::BodyId center, lupnt::BodyId target) -> lupnt::Vec6d {
            return lupnt::spice::GetBodyPosVel(t_tdb, center, target).cast<double>();
          },
          py::arg("t_tdb"), py::arg("center"), py::arg("target"))
      .def(
          "get_body_pos_vel",
          [](lupnt::VecXd t_tdb, lupnt::BodyId center, lupnt::BodyId target) -> lupnt::MatX6d {
            return lupnt::spice::GetBodyPosVel(t_tdb, center, target).cast<double>();
          },
          py::arg("t_tdb"), py::arg("center"), py::arg("target"))
      .def(
          "get_body_pos_spice",
          [](double t_tdb, lupnt::BodyId obs, lupnt::BodyId target, const std::string& ref_frame,
             const std::string& ab_orrection) -> lupnt::Vec3d {
            return lupnt::spice::GetBodyPosSpice(t_tdb, obs, target, ref_frame, ab_orrection)
                .cast<double>();
          },
          py::arg("t_tdb"), py::arg("obs"), py::arg("target"), py::arg("ref_frame"),
          py::arg("ab_correction"));
}
