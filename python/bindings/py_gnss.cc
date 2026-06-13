/**
 * @file py_gnss.cc
 * @author Stanford NAV LAB
 * @brief Python bindings for GNSS attitude / SP3 / ANTEX / RINEX-nav file
 *        interfaces.
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 * Exposes `lupnt::GnssAttitude` (`lupnt/agents/gnss_attitude.h`),
 * `lupnt::GnssYawSteering` (`lupnt/agents/gnss_yaw_steering.h` -- the GPS/
 * Galileo/BDS-3 yaw-attitude steering laws of Cheng et al. (2025),
 * https://doi.org/10.1016/j.asr.2024.10.064),
 * `lupnt::Sp3Loader` (`lupnt/interfaces/sp3_loader.h`),
 * `lupnt::AntexLoader` (`lupnt/interfaces/antex_loader.h`), and
 * `lupnt::RinexNavLoader` (`lupnt/interfaces/rinex_nav_loader.h`) to Python.
 * These are the C++ counterparts of `pylupnt.interfaces.gnss_file_loader`
 * (`SP3Loader`/`BRDCLoader`) and `pylupnt.interfaces.antex_file_loader`
 * (`ANTEXLoader`); together they implement the SP3/ANTEX phase-center-offset
 * (PCO) correction workflow demonstrated in
 * `projects/Plasmasphere_Delay_Datagen/generate_antenna_pco.ipynb` and used
 * by `GnssConstellation::SetupSatelliteStatesFromFiles`.
 *
 * Also registers the `GnssConst` / `GnssFreq` enums (`lupnt/devices/space_comms.h`)
 * needed as parameter/return types by `AntexLoader::GetPco`.
 *
 * Paths are accepted as plain strings (and converted to `std::filesystem::path`
 * internally) for consistency with `InitFile` (`py_file.cc`).
 */
#include <lupnt/lupnt.h>

#include <filesystem>
#include <string>
#include <vector>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

namespace {
  std::vector<std::filesystem::path> ToPaths(const std::vector<std::string>& filepaths) {
    std::vector<std::filesystem::path> paths;
    paths.reserve(filepaths.size());
    for (const auto& f : filepaths) paths.emplace_back(f);
    return paths;
  }
}  // namespace

void InitGnss(py::module& m) {
  // ---- Enums (shared with `AntexLoader`/`GnssConstellation`) --------------

  py::enum_<GnssConst>(m, "GnssConst")
      .value("GPS", GnssConst::GPS)
      .value("GLONASS", GnssConst::GLONASS)
      .value("GALILEO", GnssConst::GALILEO)
      .value("BEIDOU", GnssConst::BEIDOU)
      .value("QZSS", GnssConst::QZSS)
      .export_values();

  py::enum_<GnssFreq>(m, "GnssFreq")
      .value("L1", GnssFreq::L1)
      .value("L2", GnssFreq::L2)
      .value("L5", GnssFreq::L5)
      .value("E1", GnssFreq::E1)
      .value("E6", GnssFreq::E6)
      .value("E5", GnssFreq::E5)
      .value("E5a", GnssFreq::E5a)
      .value("E5b", GnssFreq::E5b)
      .export_values();

  // ---- GnssAttitude --------------------------------------------------------

  py::class_<GnssAttitude>(m, "GnssAttitude")
      .def(py::init<>())
      .def(py::init<const Vec3&, const Vec3&>(), py::arg("r_sat_eci"), py::arg("r_sun_eci"),
           "Construct and immediately compute the attitude triad (ex, ey, ez)")
      .def_static(
          "compute",
          [](const Vec3& r_sat_eci, const Vec3& r_sun_eci) {
            Vec3 ex, ey, ez;
            GnssAttitude::Compute(r_sat_eci, r_sun_eci, ex, ey, ez);
            return py::make_tuple(ex, ey, ez);
          },
          py::arg("r_sat_eci"), py::arg("r_sun_eci"),
          "Compute the orthonormal GNSS attitude triad (ex, ey, ez) directly from the "
          "Sun direction; returns (ex, ey, ez)")
      .def_static(
          "compute",
          [](const Vec3& r_sat_eci, const Vec3& v_sat_eci, const Vec3& r_sun_eci) {
            Vec3 ex, ey, ez;
            GnssAttitude::Compute(r_sat_eci, v_sat_eci, r_sun_eci, ex, ey, ez);
            return py::make_tuple(ex, ey, ez);
          },
          py::arg("r_sat_eci"), py::arg("v_sat_eci"), py::arg("r_sun_eci"),
          "Compute the orthonormal GNSS attitude triad (ex, ey, ez) via the documented "
          "nominal yaw-steering law (GnssYawSteering.nominal_yaw_angle); numerically "
          "identical to compute(r_sat_eci, r_sun_eci) but routes explicitly through the "
          "implemented yaw-steering law. Returns (ex, ey, ez)")
      .def_static(
          "compute_from_yaw_angle",
          [](const Vec3& r_sat_eci, const Vec3& v_sat_eci, Real yaw_angle) {
            Vec3 ex, ey, ez;
            GnssAttitude::ComputeFromYawAngle(r_sat_eci, v_sat_eci, yaw_angle, ex, ey, ez);
            return py::make_tuple(ex, ey, ez);
          },
          py::arg("r_sat_eci"), py::arg("v_sat_eci"), py::arg("yaw_angle"),
          "Build the attitude triad (ex, ey, ez) from an explicit yaw angle `phi` "
          "[rad] (e.g. from any GnssYawSteering nominal or modeled/maneuver yaw law) "
          "and orbital geometry alone, by rotating the orbital reference frame "
          "(ex_ref = normalize(v_sat), ey_ref = cross(ez, ex_ref)) about the nadir "
          "axis ez by `phi`. Returns (ex, ey, ez)")
      .def(
          "update",
          [](GnssAttitude& att, const Vec3& r_sat_eci, const Vec3& r_sun_eci) {
            att.Update(r_sat_eci, r_sun_eci);
          },
          py::arg("r_sat_eci"), py::arg("r_sun_eci"),
          "(Re-)compute and cache the attitude triad for this instance, directly from "
          "the Sun direction")
      .def(
          "update",
          [](GnssAttitude& att, const Vec3& r_sat_eci, const Vec3& v_sat_eci,
             const Vec3& r_sun_eci) { att.Update(r_sat_eci, v_sat_eci, r_sun_eci); },
          py::arg("r_sat_eci"), py::arg("v_sat_eci"), py::arg("r_sun_eci"),
          "(Re-)compute and cache the attitude triad for this instance via the "
          "documented nominal yaw-steering law; see compute(r_sat_eci, v_sat_eci, "
          "r_sun_eci)")
      .def("get_ex", &GnssAttitude::GetEx)
      .def("get_ey", &GnssAttitude::GetEy)
      .def("get_ez", &GnssAttitude::GetEz)
      .def("get_rotation_matrix", &GnssAttitude::GetRotationMatrix,
           "Body-to-inertial rotation matrix [ex, ey, ez] (columns)")
      .def(
          "get_angles",
          [](const GnssAttitude& att, const Vec3& u) {
            Real theta, phi;
            att.GetAngles(u, theta, phi);
            return py::make_tuple(theta, phi);
          },
          py::arg("u"),
          "Off-boresight angles (theta, phi) of unit direction `u`; returns (theta, phi)");

  // ---- GnssYawSteering ------------------------------------------------------
  // Stateless yaw-attitude steering laws (Cheng et al., 2025,
  // https://doi.org/10.1016/j.asr.2024.10.064); see `lupnt/agents/gnss_yaw_steering.h`
  // for the equation each static method mirrors. Exposed as `staticmethod`s on
  // a non-instantiable class, mirroring `AntexLoader`'s static helpers.

  py::class_<GnssYawSteering>(m, "GnssYawSteering")
      .def_static("beta_angle", &GnssYawSteering::BetaAngle, py::arg("r_sat"), py::arg("v_sat"),
                  py::arg("r_sun"), "Sun elevation angle beta above the orbital plane [rad]")
      .def_static("orbit_angle", &GnssYawSteering::OrbitAngle, py::arg("r_sat"), py::arg("v_sat"),
                  py::arg("r_sun"),
                  "Orbit angle mu [rad] from the orbit midnight point to the satellite")
      .def_static("orbit_noon_angle", &GnssYawSteering::OrbitNoonAngle, py::arg("mu"),
                  "Convert a midnight-referenced orbit angle `mu` to the noon-referenced "
                  "orbit angle `eta` (eta = wrap_to_pi(mu - pi)), see Eq. (9)")
      .def_static("nominal_yaw_angle", &GnssYawSteering::NominalYawAngle, py::arg("beta"),
                  py::arg("mu"), "Nominal yaw angle phi = atan2(-tan(beta), sin(mu))  [Eq. 1]")
      .def_static("nominal_yaw_rate", &GnssYawSteering::NominalYawRate, py::arg("beta"),
                  py::arg("mu"), py::arg("mu_dot"), "Nominal yaw rate  [Eq. 2]")
      .def_static("gps_iif_shadow_yaw_angle", &GnssYawSteering::GpsIIFShadowYawAngle, py::arg("t"),
                  py::arg("ts"), py::arg("te"), py::arg("phi_ts"), py::arg("phi_te"),
                  "GPS Block IIF constant-yaw-rate shadow-crossing model  [Eq. 3-4]")
      .def_static("gps_iif_noon_turn_yaw_angle", &GnssYawSteering::GpsIIFNoonTurnYawAngle,
                  py::arg("t"), py::arg("ts"), py::arg("phi_ts"), py::arg("beta"),
                  "GPS Block IIF noon-turn maneuver model  [Eq. 5]")
      .def_static("gps_iir_midnight_turn_yaw_angle", &GnssYawSteering::GpsIIRMidnightTurnYawAngle,
                  py::arg("t"), py::arg("ts"), py::arg("phi_ts"), py::arg("beta"),
                  "GPS Block IIR midnight-turn maneuver model  [Eq. 6]")
      .def_static("gps_iir_noon_turn_yaw_angle", &GnssYawSteering::GpsIIRNoonTurnYawAngle,
                  py::arg("t"), py::arg("ts"), py::arg("phi_ts"), py::arg("beta"),
                  "GPS Block IIR noon-turn maneuver model  [Eq. 7]")
      .def_static("gps3_sun_vector", &GnssYawSteering::Gps3SunVector, py::arg("beta"),
                  py::arg("mu"),
                  "Sun direction unit vector in the orbital reference frame, GPS III "
                  "(along-track/orbit-normal/Earth-direction) convention  [Eq. 11]")
      .def_static("gps3_eclipse_yaw_angle", &GnssYawSteering::Gps3EclipseYawAngle, py::arg("beta"),
                  py::arg("mu"), py::arg("phi_nom"),
                  "Improved rate-limited yaw-steering law for GPS Block III satellites "
                  "during eclipse-season maneuvers (Montenbruck et al., 2026)  "
                  "[Eq. 11-13, 17-18]")
      .def_static("galileo_sun_vector", &GnssYawSteering::GalileoSunVector, py::arg("eta"),
                  py::arg("beta"), "Sun reference vector in the orbital reference frame  [Eq. 9]")
      .def_static("galileo_iov_nominal_yaw_angle", &GnssYawSteering::GalileoIovNominalYawAngle,
                  py::arg("eta"), py::arg("beta"), "Galileo IOV nominal yaw angle  [Eq. 8]")
      .def_static("galileo_iov_eclipse_yaw_angle", &GnssYawSteering::GalileoIovEclipseYawAngle,
                  py::arg("eta"), py::arg("beta"), py::arg("phi_nom"),
                  "Galileo IOV modeled yaw angle during eclipse-season maneuvers  [Eq. 10-11]")
      .def_static("galileo_foc_yaw_angle", &GnssYawSteering::GalileoFocYawAngle, py::arg("t"),
                  py::arg("ts"), py::arg("phi_s"),
                  "Galileo FOC yaw-steering law during yaw maneuvers  [Eq. 12]")
      .def_static("bds3_cast_igso_yaw_angle", &GnssYawSteering::Bds3CastIgsoYawAngle, py::arg("t"),
                  py::arg("ts"), py::arg("phi_s"),
                  "BDS-3 CAST IGSO WHU yaw-steering model  [Eq. 13]")
      .def_static("bds3_cast_meo_yaw_angle", &GnssYawSteering::Bds3CastMeoYawAngle, py::arg("t"),
                  py::arg("ts"), py::arg("phi_s"),
                  "BDS-3 CAST MEO WHU yaw-steering model  [Eq. 14]")
      .def_static("bds3_secm_csno_yaw_angle", &GnssYawSteering::Bds3SecmCsnoYawAngle,
                  py::arg("beta"), py::arg("mu"),
                  "BDS-3 SECM MEO CSNO model near zero-beta yaw maneuvers  [Eq. 15]")
      .def_static("bds3_secm_mcsno_yaw_angle", &GnssYawSteering::Bds3SecmMcsnoYawAngle,
                  py::arg("t"), py::arg("t0"), py::arg("ts"), py::arg("te"), py::arg("beta"),
                  py::arg("beta_dot"), py::arg("mu"), py::arg("mu_ts"), py::arg("mu_dot"),
                  py::arg("phi_ts"),
                  "BDS-3 SECM MEO modified-CSNO (MCSNO) model with smooth zero-beta "
                  "transition  [Eq. 16]")
      .def_static("sign", &GnssYawSteering::Sign, py::arg("a"), py::arg("b"),
                  "FORTRAN SIGN(a, b) intrinsic: magnitude of `a` with the sign of `b`");

  // ---- Sp3Loader -----------------------------------------------------------

  py::class_<Sp3Loader>(m, "Sp3Loader")
      .def(py::init<>())
      .def(py::init([](const std::string& filepath) {
             return Sp3Loader(std::filesystem::path(filepath));
           }),
           py::arg("filepath"), "Construct and load a single SP3 file")
      .def(py::init([](const std::vector<std::string>& filepaths) {
             return Sp3Loader(ToPaths(filepaths));
           }),
           py::arg("filepaths"), "Construct and load multiple SP3 files (e.g. consecutive days)")
      .def(
          "load_file",
          [](Sp3Loader& loader, const std::string& filepath) {
            loader.LoadFile(std::filesystem::path(filepath));
          },
          py::arg("filepath"), "Parse an additional SP3 file and merge its samples in")
      .def("get_satellites", &Sp3Loader::GetSatellites,
           "SP3 identifiers of all satellites with data loaded so far, e.g. ['G01', 'E11', ...]")
      .def("has_satellite", &Sp3Loader::HasSatellite, py::arg("sat_id"))
      .def(
          "get_pos_vel_clock",
          [](const Sp3Loader& loader, const std::string& sat_id, Real t_tai) {
            Vec6 rv_ecef;
            Real clock_bias_s;
            loader.GetPosVelClock(sat_id, t_tai, rv_ecef, clock_bias_s);
            return py::make_tuple(rv_ecef, clock_bias_s);
          },
          py::arg("sat_id"), py::arg("t_tai"),
          "Interpolated ECEF (position, velocity) [m, m/s] and clock bias [s] "
          "(including the relativistic correction); returns (rv_ecef, clock_bias_s)")
      .def("get_pos_vel", &Sp3Loader::GetPosVel, py::arg("sat_id"), py::arg("t_tai"),
           "Interpolated ECEF position/velocity [m, m/s] only")
      .def("get_time_span", &Sp3Loader::GetTimeSpan, py::arg("sat_id"),
           "Time span [TAI seconds] covered by the loaded ephemeris for `sat_id`: (t_min, t_max)");

  // ---- AntexLoader ---------------------------------------------------------

  py::class_<AntexLoader>(m, "AntexLoader")
      .def(py::init<>())
      .def(py::init([](const std::string& filepath) {
             return AntexLoader(std::filesystem::path(filepath));
           }),
           py::arg("filepath"), "Construct and load a single ANTEX file")
      .def(
          "load_file",
          [](AntexLoader& loader, const std::string& filepath) {
            loader.LoadFile(std::filesystem::path(filepath));
          },
          py::arg("filepath"), "Parse an additional ANTEX file and merge its satellite entries in")
      .def("get_pco",
           py::overload_cast<GnssConst, int, GnssFreq, Real>(&AntexLoader::GetPco, py::const_),
           py::arg("gnss_const"), py::arg("prn"), py::arg("freq"), py::arg("t_tai"),
           "Satellite antenna phase-center offset (PCO), in the satellite North/East/Up "
           "(NEU) antenna frame [m]")
      .def("get_pco",
           py::overload_cast<const std::string&, int, const std::string&, Real>(
               &AntexLoader::GetPco, py::const_),
           py::arg("gnss_letter"), py::arg("prn"), py::arg("antex_freq_code"), py::arg("t_tai"),
           "Same as `get_pco`, but specifying the raw ANTEX frequency code directly "
           "(e.g. 'G01', 'E05')")
      .def("get_available_freq_codes", &AntexLoader::GetAvailableFreqCodes, py::arg("gnss_const"),
           py::arg("prn"), py::arg("t_tai"),
           "ANTEX frequency codes available for satellite `prn` at epoch `t_tai`")
      .def("has_satellite", &AntexLoader::HasSatellite, py::arg("gnss_const"), py::arg("prn"))
      .def_static(
          "compute_ijk_to_ecef_rotation", &AntexLoader::ComputeIjkToEcefRotation, py::arg("t_tai"),
          py::arg("r_sat_ecef"),
          "'IJK' rotation matrix Cijk s.t. pco_ecef = Cijk @ pco_neu "
          "(ivec = jvec x kvec, jvec = normalize(r_sun_ecef - r_sat_ecef), "
          "kvec = -normalize(r_sat_ecef); intentionally distinct from "
          "GnssAttitude's orthonormal body frame, see `phase_center_offset.py::ijk_to_ecef_rot`)")
      .def_static("apply_pco_correction_ecef", &AntexLoader::ApplyPcoCorrectionEcef,
                  py::arg("t_tai"), py::arg("pos_sp3_ecef"), py::arg("pco_neu_m"),
                  "pos_corrected = pos_sp3_ecef + Cijk(t_tai, pos_sp3_ecef) @ pco_neu "
                  "(mirrors the workflow in generate_antenna_pco.ipynb)")
      .def_static("gnss_letter", &AntexLoader::GnssLetter, py::arg("gnss_const"),
                  "Single-letter GNSS system identifier, e.g. GnssConst.GPS -> 'G'")
      .def_static("sat_id", &AntexLoader::SatId, py::arg("gnss_const"), py::arg("prn"),
                  "SP3/ANTEX/RINEX satellite identifier, e.g. (GnssConst.GPS, 5) -> 'G05'");

  // ---- RinexNavLoader ------------------------------------------------------

  py::class_<RinexNavLoader>(m, "RinexNavLoader")
      .def(py::init<>())
      .def(py::init([](const std::string& filepath) {
             return RinexNavLoader(std::filesystem::path(filepath));
           }),
           py::arg("filepath"), "Construct and load a single RINEX nav file")
      .def(py::init([](const std::vector<std::string>& filepaths) {
             return RinexNavLoader(ToPaths(filepaths));
           }),
           py::arg("filepaths"),
           "Construct and load multiple RINEX nav files (e.g. consecutive days)")
      .def(
          "load_file",
          [](RinexNavLoader& loader, const std::string& filepath) {
            loader.LoadFile(std::filesystem::path(filepath));
          },
          py::arg("filepath"), "Parse an additional RINEX nav file and merge its messages in")
      .def("get_satellites", &RinexNavLoader::GetSatellites,
           "Identifiers of all satellites with navigation messages loaded so far "
           "(GLONASS excluded), e.g. ['G01', 'E11', ...]")
      .def("has_satellite", &RinexNavLoader::HasSatellite, py::arg("sat_id"))
      .def(
          "get_pos_vel_clock",
          [](const RinexNavLoader& loader, const std::string& sat_id, Real t_tai) {
            Vec6 rv_ecef;
            Real clock_corr_s;
            loader.GetPosVelClock(sat_id, t_tai, rv_ecef, clock_corr_s);
            return py::make_tuple(rv_ecef, clock_corr_s);
          },
          py::arg("sat_id"), py::arg("t_tai"),
          "Broadcast-ephemeris ECEF (position, velocity) [m, m/s] and clock correction [s] "
          "via Keplerian propagation of the closest navigation message; "
          "returns (rv_ecef, clock_corr_s)")
      .def("get_pos_vel", &RinexNavLoader::GetPosVel, py::arg("sat_id"), py::arg("t_tai"),
           "Broadcast-ephemeris ECEF position/velocity [m, m/s] only");
}
