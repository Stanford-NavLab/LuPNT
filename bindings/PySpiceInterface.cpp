#include <lupnt/physics/SpiceInterface.h>
#include <pybind11/pybind11.h>
<<<<<<< HEAD

#include "Matrix.hpp"
=======
>>>>>>> 0c6a8c8f4d4a0d22d81dfc2122faa380efd77274

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module &m) {
  py::class_<SpiceInterface>(m, "SpiceInterface")
      .def_static("load_spice_kernel", &LPT::SpiceInterface::LoadSpiceKernel)
      .def_static("extract_pck_coeffs", &LPT::SpiceInterface::ExtractPckCoeffs)
<<<<<<< HEAD
=======
      .def_static("get_frame_conversion_matrix",
                  [](double et, std::string from, std::string to) {
                    return LPT::SpiceInterface::GetFrameConversionMatrix(
                        et, from, to);
                  })

>>>>>>> 0c6a8c8f4d4a0d22d81dfc2122faa380efd77274
      .def_static("string_to_tdb", py::overload_cast<std::string>(
                                       &LPT::SpiceInterface::StringToTDB))
      .def_static("string_to_tai", py::overload_cast<std::string>(
                                       &LPT::SpiceInterface::StringToTAI))
      .def_static("tdb_to_string_utc",
<<<<<<< HEAD
                  py::overload_cast<ad::real, int>(
                      &LPT::SpiceInterface::TDBtoStringUTC))
      .def_static("convert_time",
                  py::overload_cast<ad::real, std::string, std::string>(
                      &LPT::SpiceInterface::ConvertTime))
      .def_static("get_body_pos_vel", py::overload_cast<ad::real, int, int>(
                                          &LPT::SpiceInterface::GetBodyPosVel))
      .def_static(
          "get_body_pos",
          py::overload_cast<std::string, ad::real, std::string, std::string,
                            std::string>(&LPT::SpiceInterface::GetBodyPos))
      .def_static(
          "get_frame_conversion_matrix",
          [](ad::real et, std::string from_frame, std::string to_frame) {
            const Eigen::MatrixXd M_rot =
                LPT::SpiceInterface::GetFrameConversionMatrix(et, from_frame,
                                                              to_frame);
            return EigenMatrixToPythonArray(M_rot);
          });
=======
                  [](double tdb, int prec) {
                    return LPT::SpiceInterface::TDBtoStringUTC(tdb, prec);
                  })
      .def_static("convert_time",
                  [](double et, std::string from, std::string to) {
                    return LPT::SpiceInterface::ConvertTime(et, from, to);
                  })
      .def_static("get_body_pos_vel",
                  [](double ta, int center, int target) {
                    return LPT::SpiceInterface::GetBodyPosVel(ta, center,
                                                              target);
                  })
      .def_static("get_body_pos", [](std::string targetName, ad::real epoch,
                                     std::string refFrame, std::string obsName,
                                     std::string abCorrection) {
        return LPT::SpiceInterface::GetBodyPos(targetName, epoch, refFrame,
                                               obsName, abCorrection);
      });
>>>>>>> 0c6a8c8f4d4a0d22d81dfc2122faa380efd77274
}