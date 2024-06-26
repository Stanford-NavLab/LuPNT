#include <lupnt/core/constants.h>
#include <lupnt/physics/spice_interface.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module &m) {
  py::class_<SpiceInterface>(m, "SpiceInterface")
      .def_static("load_spice_kernel", &lupnt::LoadSpiceKernel)
      .def_static("extract_pck_coeffs", &lupnt::ExtractPckCoeffs)
      .def_static("get_frame_conversion_matrix",
                  [](double et, std::string from, std::string to) {
                    return lupnt::GetFrameConversionMat(et, from, to);
                  })

      .def_static(
          "string_to_tdb",
          [](std::string utc) -> double {
            return lupnt::StringToTDB(utc).val();
          },
          py::arg("gregorian_date"))
      .def_static(
          "string_to_tai",
          [](std::string utc) -> double {
            return lupnt::StringToTAI(utc).val();
          },
          py::arg("gregorian_date"))
      .def_static(
          "tdb_to_string_utc",
          [](double tdb, int prec) { return lupnt::TDBtoStringUTC(tdb, prec); })
      .def_static("convert_time",
                  [](double et, std::string from, std::string to) {
                    return lupnt::ConvertTime(et, from, to);
                  })
      .def_static("get_body_pos_vel",
                  [](double ta, int center, int target) {
                    return lupnt::GetBodyPosVel(ta, center, target);
                  })
      .def_static("get_body_pos", [](std::string targetName, lupnt::real epoch,
                                     std::string refFrame, std::string obsName,
                                     std::string abCorrection) {
        return lupnt::GetBodyPos(targetName, epoch, refFrame, obsName,
                                 abCorrection);
      });
}