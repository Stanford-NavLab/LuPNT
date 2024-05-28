#include <lupnt/core/constants.h>
#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/spice_interface.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module &m) {
  py::class_<SpiceInterface>(m, "SpiceInterface")
      .def_static("load_spice_kernel", &lupnt::SpiceInterface::LoadSpiceKernel)
      .def_static("extract_pck_coeffs",
                  &lupnt::SpiceInterface::ExtractPckCoeffs)
      .def_static(
          "get_frame_conversion_matrix",
          [](double t_tai, lupnt::Frame from,
             lupnt::Frame to) -> lupnt::Matrix6d {
            return lupnt::SpiceInterface::GetFrameConversionMatrix(t_tai, from,
                                                                   to)
                .cast<double>();
          },
          py::arg("t_tai"), py::arg("from"), py::arg("to"))
      .def_static(
          "string_to_tdb",
          [](std::string utc) -> double {
            return lupnt::SpiceInterface::StringToTDB(utc).val();
          },
          py::arg("gregorian_date"))
      .def_static(
          "string_to_tai",
          [](std::string utc) -> double {
            return lupnt::SpiceInterface::StringToTAI(utc).val();
          },
          py::arg("gregorian_date"))
      .def_static(
          "tdb_to_string_utc",
          [](double t_tdb, int prec) {
            return lupnt::SpiceInterface::TDBtoStringUTC(t_tdb, prec);
          },
          py::arg("t_tdb"), py::arg("precision"))
      .def_static(
          "convert_time",
          [](double t_tai, std::string from, std::string to) {
            return lupnt::SpiceInterface::ConvertTime(t_tai, from, to).val();
          },
          py::arg("t_tai"), py::arg("from"), py::arg("to"))
      .def_static(
          "get_body_pos_vel",
          [](double t_tai, lupnt::NaifId center,
             lupnt::NaifId target) -> lupnt::Vector6d {
            return lupnt::SpiceInterface::GetBodyPosVel(t_tai, center, target)
                .cast<double>();
          },
          py::arg("t_tai"), py::arg("center"), py::arg("target"))
      .def_static(
          "get_body_pos_vel",
          [](lupnt::VectorXd t_tai, lupnt::NaifId center,
             lupnt::NaifId target) -> lupnt::Matrixd<-1, 6> {
            return lupnt::SpiceInterface::GetBodyPosVel(t_tai, center, target)
                .cast<double>();
          },
          py::arg("t_tai"), py::arg("center"), py::arg("target"))
      .def_static(
          "get_body_pos",
          [](lupnt::NaifId targetName, double t_tai, lupnt::Frame refFrame,
             lupnt::NaifId obsName,
             std::string abCorrection) -> lupnt::Vector3d {
            return lupnt::SpiceInterface::GetBodyPos(
                       targetName, t_tai, refFrame, obsName, abCorrection)
                .cast<double>();
          },
          py::arg("target_name"), py::arg("t_tai"), py::arg("ref_frame"),
          py::arg("obs_name"), py::arg("ab_correction"));
}