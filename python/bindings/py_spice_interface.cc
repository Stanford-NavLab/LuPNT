#include <lupnt/lupnt.h>

// Enables the std::filesystem::path <-> str caster so load_spice_kernel() can
// be called with a plain Python string.
#include <pybind11/stl/filesystem.h>

#include "py_pybind11.h"

namespace py = pybind11;

class SpiceInterface {};  // dummy class

void init_spice_interface(py::module& m) {
  auto m_spice = m.def_submodule("spice", "This is A.");

  m_spice
      .def("load_spice_kernel",
           static_cast<void (*)(const std::filesystem::path&)>(&lupnt::spice::LoadSpiceKernel),
           "Load an additional SPICE kernel file (loading LuPNT's default kernel set first if "
           "needed).")
      .def("extract_pck_coeffs", &lupnt::spice::ExtractPckCoeffs,
           "Debug helper: extract and print lunar PCK orientation (Euler-angle) coefficients.")
      .def(
          "get_frame_conversion_mat",
          [](double t_tdb, const std::string& from, const std::string& to) -> lupnt::Mat6d {
            return lupnt::spice::GetFrameConversionMat(t_tdb, from, to).cast<double>();
          },
          py::arg("t_tdb"), py::arg("from"), py::arg("to"),
          "6x6 state (position+velocity) rotation matrix from SPICE frame `from` to `to` at "
          "`t_tdb` [s, TDB].")
      .def(
          "string2tdb",
          [](std::string utc) -> double { return lupnt::spice::StringToTdb(utc).val(); },
          py::arg("gregorian_date"),
          "Parse a calendar/Julian-date string to TDB seconds past J2000.")
      .def(
          "string2tai",
          [](std::string utc) -> double { return lupnt::spice::StringToTai(utc).val(); },
          py::arg("gregorian_date"),
          "Parse a calendar/Julian-date string to TAI seconds past J2000.")
      .def(
          "tdb2string_utc",
          [](double t_tdb, int prec) { return lupnt::spice::TDBtoStringUTC(t_tdb, prec); },
          py::arg("t_tdb"), py::arg("precision"),
          "Format a TDB epoch [s past J2000] as a UTC calendar string (`precision` "
          "fractional-second digits).")
      .def("convert_time", &lupnt::spice::ConvertTime, py::arg("t"), py::arg("from"), py::arg("to"),
           "Convert a time value between time systems (TAI, TDB, TT, UTC, GPS, ...) via SPICE [s].")
      .def(
          "get_time_ephemeris_offset",
          [](double t_tdb, int naif_id) {
            return lupnt::spice::GetTimeEphemerisOffset(t_tdb, naif_id).val();
          },
          py::arg("t_tdb"), py::arg("naif_id"),
          "Scalar time-ephemeris offset [s] at TDB epoch `t_tdb` from a loaded SPK time-ephemeris "
          "segment. Known ids: 1000000001 = TT-TDB (de440t.bsp), 1000000005 = TDB-TCL, "
          "1000000006 = TDB-TL (lunar time ephemerides).")
      .def(
          "get_time_ephemeris_offset",
          [](const lupnt::VecXd& t_tdb, int naif_id) {
            lupnt::VecXd out(t_tdb.size());
            for (int i = 0; i < t_tdb.size(); i++)
              out(i) = lupnt::spice::GetTimeEphemerisOffset(t_tdb(i), naif_id).val();
            return out;
          },
          py::arg("t_tdb"), py::arg("naif_id"),
          "Vectorized scalar time-ephemeris offset [s] at each TDB epoch.")
      .def(
          "get_body_pos_vel",
          [](double t_tdb, lupnt::BodyId center, lupnt::BodyId target) -> lupnt::Vec6d {
            return lupnt::spice::GetBodyPosVel(t_tdb, center, target).cast<double>();
          },
          py::arg("t_tdb"), py::arg("center"), py::arg("target"),
          "Inertial (J2000) position+velocity of `target` relative to `center` at `t_tdb` "
          "[km, km/s].")
      .def(
          "get_body_pos_vel",
          [](lupnt::VecXd t_tdb, lupnt::BodyId center, lupnt::BodyId target) -> lupnt::MatX6d {
            return lupnt::spice::GetBodyPosVel(t_tdb, center, target).cast<double>();
          },
          py::arg("t_tdb"), py::arg("center"), py::arg("target"),
          "Vectorized: per-epoch J2000 position+velocity of `target` relative to `center` "
          "[N x 6, km & km/s].")
      .def(
          "get_body_pos_spice",
          [](double t_tdb, lupnt::BodyId obs, lupnt::BodyId target, const std::string& ref_frame,
             const std::string& ab_orrection) -> lupnt::Vec3d {
            return lupnt::spice::GetBodyPosSpice(t_tdb, obs, target, ref_frame, ab_orrection)
                .cast<double>();
          },
          py::arg("t_tdb"), py::arg("obs"), py::arg("target"), py::arg("ref_frame"),
          py::arg("ab_correction"),
          "Position of `target` relative to `obs` at `t_tdb` from SPICE (spkpos_c), in "
          "`ref_frame` with aberration correction `ab_correction` [km].")
      .def(
          "get_body_pos_vel_spice",
          [](double t_tdb, lupnt::BodyId obs, lupnt::BodyId target, const std::string& ref_frame,
             const std::string& ab_correction) -> lupnt::Vec6d {
            return lupnt::spice::GetBodyPosVelSpice(t_tdb, obs, target, ref_frame, ab_correction);
          },
          py::arg("t_tdb"), py::arg("obs"), py::arg("target"), py::arg("ref_frame") = "J2000",
          py::arg("ab_correction") = "NONE",
          "Position+velocity of `target` relative to `obs` at `t_tdb` from SPICE (spkezr_c), in "
          "`ref_frame` with aberration correction `ab_correction` [km, km/s]. Unlike "
          "`get_body_pos_vel` this reads the furnished SPK kernels directly, so it works for "
          "spacecraft bodies loaded via `load_spice_kernel`.");
}
