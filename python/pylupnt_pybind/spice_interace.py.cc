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
          [](double et, lupnt::Frame from, lupnt::Frame to) -> lupnt::Matrix6d {
            return lupnt::SpiceInterface::GetFrameConversionMatrix(et, from, to)
                .cast<double>();
          })

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
      .def_static("tdb_to_string_utc",
                  [](double tdb, int prec) {
                    return lupnt::SpiceInterface::TDBtoStringUTC(tdb, prec);
                  })
      .def_static(
          "convert_time",
          [](double et, std::string from, std::string to) {
            return lupnt::SpiceInterface::ConvertTime(et, from, to).val();
          })
      .def_static("get_body_pos_vel",
                  [](double ta, lupnt::NaifId center,
                     lupnt::NaifId target) -> lupnt::Vector6d {
                    return lupnt::SpiceInterface::GetBodyPosVel(ta, center,
                                                                target)
                        .cast<double>();
                  })
      .def_static("get_body_pos_vel",
                  [](lupnt::VectorXd tai, lupnt::NaifId center,
                     lupnt::NaifId target) -> lupnt::Matrixd<-1, 6> {
                    return lupnt::SpiceInterface::GetBodyPosVel(tai, center,
                                                                target)
                        .cast<double>();
                  })
      .def_static("get_body_pos",
                  [](lupnt::NaifId targetName, lupnt::real epoch,
                     lupnt::Frame refFrame, lupnt::NaifId obsName,
                     std::string abCorrection) -> lupnt::Vector3d {
                    return lupnt::SpiceInterface::GetBodyPos(targetName, epoch,
                                                             refFrame, obsName,
                                                             abCorrection)
                        .cast<double>();
                  });
}