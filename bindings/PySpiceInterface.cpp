#include <pybind11/pybind11.h>
#include <lupnt/physics/SpiceInterface.h>

namespace py = pybind11;

class SpiceInterface {}; //dummy class

void init_spice_interface(py::module &m)
{
  py::class_<SpiceInterface>(m, "SpiceInterface")
    .def_static("load_spice_kernel", &LPT::SpiceInterface::LoadSpiceKernel)
    .def_static("extract_pck_coeffs", &LPT::SpiceInterface::ExtractPckCoeffs)
    .def_static("get_frame_conversion_matrix", py::overload_cast<ad::real, std::string, std::string>(&LPT::SpiceInterface::GetFrameConversionMatrix))
    .def_static("string_to_tdb", py::overload_cast<std::string>(&LPT::SpiceInterface::StringToTDB))
    .def_static("string_to_tai", py::overload_cast<std::string>(&LPT::SpiceInterface::StringToTAI))
    .def_static("tdb_to_string_utc", py::overload_cast<ad::real, int>(&LPT::SpiceInterface::TDBtoStringUTC))
    .def_static("convert_time", py::overload_cast<ad::real, std::string, std::string>(&LPT::SpiceInterface::ConvertTime))
    .def_static("get_body_pos_vel", py::overload_cast<ad::real, int, int>(&LPT::SpiceInterface::GetBodyPosVel))
    .def_static("get_body_pos", py::overload_cast<std::string, ad::real, std::string, std::string, std::string>(&LPT::SpiceInterface::GetBodyPos));
}