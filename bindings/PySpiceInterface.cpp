#include <lupnt/physics/spice_interface.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module &m) {
  py::class_<SpiceInterface>(m, "SpiceInterface")
      .def_static("load_spice_kernel", &lupnt::SpiceInterface::LoadSpiceKernel)
      .def_static("extract_pck_coeffs", &lupnt::SpiceInterface::ExtractPckCoeffs)
      .def_static("get_frame_conversion_matrix",
                  [](double et, std::string from, std::string to) {
                    return lupnt::SpiceInterface::GetFrameConversionMatrix(
                        et, from, to);
                  })

      .def_static("string_to_tdb", py::overload_cast<std::string>(
                                       &lupnt::SpiceInterface::StringToTDB))
      .def_static("string_to_tai", py::overload_cast<std::string>(
                                       &lupnt::SpiceInterface::StringToTAI))
      .def_static("tdb_to_string_utc",
                  [](double tdb, int prec) {
                    return lupnt::SpiceInterface::TDBtoStringUTC(tdb, prec);
                  })
      .def_static("convert_time",
                  [](double et, std::string from, std::string to) {
                    return lupnt::SpiceInterface::ConvertTime(et, from, to);
                  })
      .def_static("get_body_pos_vel",
                  [](double ta, int center, int target) {
                    return lupnt::SpiceInterface::GetBodyPosVel(ta, center,
                                                              target);
                  })
      .def_static("get_body_pos", [](std::string targetName, ad::real epoch,
                                     std::string refFrame, std::string obsName,
                                     std::string abCorrection) {
        return lupnt::SpiceInterface::GetBodyPos(targetName, epoch, refFrame,
                                               obsName, abCorrection);
      });
}