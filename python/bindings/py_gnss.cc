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
 * Exposes `lupnt::GnssAttitude` (`lupnt/attitude/gnss_attitude.h`),
 * `lupnt::GnssYawSteering` (`lupnt/attitude/gnss_yaw_steering.h` -- the GPS/
 * Galileo/BDS-3 yaw-attitude steering laws of Cheng et al. (2025),
 * https://doi.org/10.1016/j.asr.2024.10.064),
 * `lupnt::Sp3Loader` (`lupnt/interfaces/sp3_loader.h`),
 * `lupnt::AntexLoader` (`lupnt/interfaces/antex_loader.h`), and
 * `lupnt::RinexNavLoader` (`lupnt/interfaces/rinex_nav_loader.h`) to Python.
 * These are the C++ counterparts of the legacy Python SP3/BRDC/ANTEX loaders
 * (`SP3Loader`/`BRDCLoader`/`ANTEXLoader`, now living under
 * `projects/Plasmasphere_Delay_Datagen/src/interfaces/`); together they
 * implement the SP3/ANTEX phase-center-offset (PCO) correction workflow
 * demonstrated in `projects/Plasmasphere_Delay_Datagen/generate_antenna_pco.ipynb`
 * and used by `GnssConstellation::SetupSatelliteStatesFromFiles`.
 *
 * Also registers the `GnssConst` / `GnssFreq` enums (`lupnt/devices/gnss_device.h`)
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

  py::enum_<GnssConst>(m, "GnssConst", "GNSS constellation / system identifier")
      .value("GPS", GnssConst::GPS, "Global Positioning System (USA)")
      .value("GLONASS", GnssConst::GLONASS, "GLONASS (Russia)")
      .value("GALILEO", GnssConst::GALILEO, "Galileo (European Union)")
      .value("BEIDOU", GnssConst::BEIDOU, "BeiDou (China)")
      .value("QZSS", GnssConst::QZSS, "Quasi-Zenith Satellite System (Japan)")
      .export_values();

  py::enum_<GnssFreq>(m, "GnssFreq", "GNSS carrier frequency band identifier")
      .value("L1", GnssFreq::L1, "GPS/QZSS L1 band (1575.42 MHz)")
      .value("L2", GnssFreq::L2, "GPS/QZSS L2 band (1227.60 MHz)")
      .value("L5", GnssFreq::L5, "GPS/QZSS L5 band (1176.45 MHz)")
      .value("E1", GnssFreq::E1, "Galileo E1 band (1575.42 MHz)")
      .value("E6", GnssFreq::E6, "Galileo E6 band (1278.75 MHz)")
      .value("E5", GnssFreq::E5, "Galileo E5 (E5a+E5b) band (1191.795 MHz)")
      .value("E5a", GnssFreq::E5a, "Galileo E5a band (1176.45 MHz)")
      .value("E5b", GnssFreq::E5b, "Galileo E5b band (1207.14 MHz)")
      .export_values();

  // ---- GnssAttitude --------------------------------------------------------

  py::class_<GnssAttitude>(m, "GnssAttitude",
                           "GNSS satellite nominal (Sun-pointing yaw-steering) attitude: computes "
                           "and caches the orthonormal body triad (ex, ey, ez) from orbit and Sun "
                           "geometry")
      .def(py::init<>(), "Construct an empty attitude (triad uninitialized until compute/update)")
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
      .def("get_ex", &GnssAttitude::GetEx,
           "Cached along-track body axis ex (unit vector, inertial frame)")
      .def("get_ey", &GnssAttitude::GetEy,
           "Cached cross-track (Sun-side) body axis ey (unit vector, inertial frame)")
      .def("get_ez", &GnssAttitude::GetEz,
           "Cached nadir body axis ez (unit vector, inertial frame)")
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
  // https://doi.org/10.1016/j.asr.2024.10.064); see `lupnt/attitude/gnss_yaw_steering.h`
  // for the equation each static method mirrors. Exposed as `staticmethod`s on
  // a non-instantiable class, mirroring `AntexLoader`'s static helpers.

  py::class_<GnssYawSteering>(
      m, "GnssYawSteering",
      "Stateless GNSS yaw-attitude steering laws (Cheng et al., 2025, "
      "https://doi.org/10.1016/j.asr.2024.10.064); non-instantiable, all methods are static")
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

  py::class_<Sp3Loader>(m, "Sp3Loader",
                        "Loader for SP3 precise-ephemeris files: parses and interpolates GNSS "
                        "satellite ECEF position/velocity and clock bias")
      .def(py::init<>(), "Construct an empty loader (no files loaded)")
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
      .def("has_satellite", &Sp3Loader::HasSatellite, py::arg("sat_id"),
           "True if ephemeris data for `sat_id` (e.g. 'G01') has been loaded")
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
           "Time span [TAI seconds] covered by the loaded ephemeris for `sat_id`: (t_min, t_max)")
      .def_static(
          "download_file_for_epoch",
          [](Real epoch, Time time_scale, const std::string& cache_dir) {
            return Sp3Loader::DownloadFileForEpoch(epoch, time_scale,
                                                   cache_dir.empty()
                                                       ? std::filesystem::path()
                                                       : std::filesystem::path(cache_dir))
                .string();
          },
          py::arg("epoch"), py::arg("time_scale") = Time::UTC, py::arg("cache_dir") = std::string(),
          "Download/cache the COD MGEX final SP3 product covering `epoch` (in `time_scale` "
          "seconds) "
          "from NASA CDDIS and return the local uncompressed .SP3 path. Needs Earthdata "
          "credentials "
          "(~/.netrc or EARTHDATA_USERNAME/EARTHDATA_PASSWORD). Cache defaults to "
          "output/gnss_files/sp3.")
      .def_static("filename_for_epoch", &Sp3Loader::FilenameForEpoch, py::arg("epoch"),
                  py::arg("time_scale") = Time::UTC,
                  "SP3 filename covering `epoch`, e.g. 'COD0MGXFIN_20260140000_01D_05M_ORB.SP3'")
      .def_static("url_for_epoch", &Sp3Loader::UrlForEpoch, py::arg("epoch"),
                  py::arg("time_scale") = Time::UTC, "CDDIS download URL for the SP3 product");

  // ---- AntexLoader ---------------------------------------------------------

  py::class_<AntexLoader>(m, "AntexLoader",
                          "Loader for ANTEX antenna files: provides satellite antenna phase-center "
                          "offsets (PCO) and the SP3 PCO-correction workflow")
      .def(py::init<>(), "Construct an empty loader (no files loaded)")
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
      .def("has_satellite", &AntexLoader::HasSatellite, py::arg("gnss_const"), py::arg("prn"),
           "True if an ANTEX entry exists for satellite (`gnss_const`, `prn`)")
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

  py::class_<RinexNavLoader>(m, "RinexNavLoader",
                             "Loader for RINEX navigation (broadcast-ephemeris) files: Keplerian "
                             "propagation of GNSS satellite ECEF state and clock correction")
      .def(py::init<>(), "Construct an empty loader (no files loaded)")
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
      .def("has_satellite", &RinexNavLoader::HasSatellite, py::arg("sat_id"),
           "True if a navigation message for `sat_id` (e.g. 'G01') has been loaded")
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
           "Broadcast-ephemeris ECEF position/velocity [m, m/s] only")
      .def_static(
          "download_file_for_epoch",
          [](Real epoch, Time time_scale, const std::string& cache_dir) {
            return RinexNavLoader::DownloadFileForEpoch(epoch, time_scale,
                                                        cache_dir.empty()
                                                            ? std::filesystem::path()
                                                            : std::filesystem::path(cache_dir))
                .string();
          },
          py::arg("epoch"), py::arg("time_scale") = Time::UTC, py::arg("cache_dir") = std::string(),
          "Download/cache the BRDC00IGS broadcast-ephemeris RINEX nav file covering `epoch` "
          "(in `time_scale` seconds) from NASA CDDIS and return the local uncompressed .rnx path. "
          "Needs Earthdata credentials (~/.netrc or EARTHDATA_USERNAME/EARTHDATA_PASSWORD). "
          "Cache defaults to output/gnss_files/brdc.")
      .def_static("filename_for_epoch", &RinexNavLoader::FilenameForEpoch, py::arg("epoch"),
                  py::arg("time_scale") = Time::UTC,
                  "BRDC filename covering `epoch`, e.g. 'BRDC00IGS_R_20260140000_01D_MN.rnx'");

  // ---- GnssReceiverParams ----------------------------------------------------

  py::class_<GnssReceiverParams>(m, "GnssReceiverParams",
                                 "Receiver tracking-loop and link-budget parameters used to derive "
                                 "per-channel measurement-noise sigmas and CN0")
      .def(py::init<>(), "Construct with default receiver parameter values")
      .def_readwrite("Bp", &GnssReceiverParams::Bp, "Carrier loop noise bandwidth [Hz]")
      .def_readwrite("T", &GnssReceiverParams::T, "Tracking loop integration time [s]")
      .def_readwrite("b", &GnssReceiverParams::b, "Front-end bandwidth factor")
      .def_readwrite("Bn", &GnssReceiverParams::Bn, "Code loop noise bandwidth [Hz]")
      .def_readwrite("Bf", &GnssReceiverParams::Bf, "Frequency loop noise bandwidth [Hz]")
      .def_readwrite("D", &GnssReceiverParams::D, "Early-to-late correlator spacing [chip]")
      .def_readwrite("L_ad", &GnssReceiverParams::L_ad, "A/D converter loss [dB]")
      .def_readwrite("L_pol", &GnssReceiverParams::L_pol, "Polarization loss [dB]")
      .def_readwrite("L_atm", &GnssReceiverParams::L_atm, "Atmospheric loss [dB]")
      .def_readwrite("T_eff", &GnssReceiverParams::T_eff, "Effective noise temperature [K]");

  // ---- GnssOccludingBody -----------------------------------------------------

  py::class_<GnssOccludingBody>(m, "GnssOccludingBody",
                                "Spherical body (e.g. Earth or Moon) that can occlude the "
                                "transmitter-receiver line of sight during visibility checks")
      .def(py::init<>(), "Construct a zero-radius body at the origin")
      .def_readwrite("radius_m", &GnssOccludingBody::radius_m, "Body radius [m]")
      .def_property(
          "position_m",
          [](const GnssOccludingBody& b) -> VecXd { return b.position_m.cast<double>(); },
          [](GnssOccludingBody& b, const VecXd& v) { b.position_m = v.cast<Real>(); },
          "Body center [m] in the same frame as satellite/receiver states (used when "
          "position_provider is None)")
      .def_property(
          "position_provider",
          [](const GnssOccludingBody& b) -> py::object {
            if (!b.position_provider) return py::none();
            return py::cpp_function([fn = b.position_provider](double t) -> VecXd {
              return fn(static_cast<Real>(t)).cast<double>();
            });
          },
          [](GnssOccludingBody& b, py::object fn) {
            if (fn.is_none()) {
              b.position_provider = nullptr;
            } else {
              b.position_provider = [fn](Real t) -> Vec3 {
                py::gil_scoped_acquire gil;
                return fn(static_cast<double>(t)).cast<Vec3>();
              };
            }
          },
          "Per-epoch position callback(t_tai) -> [x,y,z] [m]; overrides position_m if set");

  // ---- GnssMeasurementOptions ------------------------------------------------

  py::class_<GnssMeasurementOptions>(
      m, "GnssMeasurementOptions",
      "Configuration for GNSS measurement construction: frames, time "
      "scales, and light-time/relativity/visibility/CN0 modeling "
      "toggles used by BuildChannels()/Compute()")
      .def(py::init<>(), "Construct with default measurement options")
      .def_readwrite("frame", &GnssMeasurementOptions::frame,
                     "Frame in which receiver and transmitter states are expressed")
      .def_readwrite("receive_time_scale", &GnssMeasurementOptions::receive_time_scale,
                     "Time scale of the receiver signal-reception epochs passed to "
                     "Compute()/Precompute()")
      .def_readwrite("ephemeris_time_scale", &GnssMeasurementOptions::ephemeris_time_scale,
                     "Time scale in which the constellation ephemeris epochs are represented")
      .def_readwrite("solve_light_time", &GnssMeasurementOptions::solve_light_time,
                     "Iteratively solve the transmit epoch for signal light-time delay")
      .def_readwrite("apply_transmitter_relativity",
                     &GnssMeasurementOptions::apply_transmitter_relativity,
                     "Apply the transmitter special-relativistic clock correction")
      .def_readwrite("apply_shapiro_delay", &GnssMeasurementOptions::apply_shapiro_delay,
                     "Apply the Shapiro (gravitational) signal-propagation delay [m]")
      .def_readwrite("apply_visibility", &GnssMeasurementOptions::apply_visibility,
                     "Drop channels occluded by the configured occluding bodies")
      .def_readwrite("apply_cn0_threshold", &GnssMeasurementOptions::apply_cn0_threshold,
                     "Drop channels whose CN0 falls below the acquisition/tracking thresholds")
      .def_readwrite("cn0_threshold_dbhz", &GnssMeasurementOptions::cn0_threshold_dbhz,
                     "Deprecated single CN0 threshold [dBHz]; use the acquisition/tracking "
                     "thresholds below")
      .def_readwrite("cn0_acquisition_threshold_dbhz",
                     &GnssMeasurementOptions::cn0_acquisition_threshold_dbhz,
                     "Min CN0 to acquire a new satellite [dBHz]")
      .def_readwrite("cn0_tracking_threshold_dbhz",
                     &GnssMeasurementOptions::cn0_tracking_threshold_dbhz,
                     "Min CN0 to maintain an existing lock [dBHz]")
      .def_readwrite("apply_ionosphere_plasma_delay",
                     &GnssMeasurementOptions::apply_ionosphere_plasma_delay,
                     "Apply the ionosphere/plasmasphere signal delay [m]");

  // ---- GnssChannel -----------------------------------------------------------

  py::class_<GnssChannel>(
      m, "GnssChannel",
      "A single transmitter-receiver GNSS signal channel: transmitter "
      "identity/state plus per-channel delays, CN0, and observable noise sigmas")
      .def(py::init<>(), "Construct an empty channel with default field values")
      .def_readwrite("gnss_const", &GnssChannel::gnss_const, "Transmitter GNSS constellation")
      .def_readwrite("prn", &GnssChannel::prn, "Transmitter satellite PRN")
      .def_readwrite("frequency", &GnssChannel::frequency, "Carrier frequency band of this channel")
      .def_readwrite("receive_time", &GnssChannel::receive_time,
                     "Signal-reception epoch [s] in receive_time_scale")
      .def_readwrite("transmit_time", &GnssChannel::transmit_time,
                     "Signal-transmission epoch [s] in transmit_time_scale")
      .def_property(
          "tx_state", [](const GnssChannel& ch) -> VecXd { return ch.tx_state.cast<double>(); },
          [](GnssChannel& ch, const VecXd& v) { ch.tx_state = v.cast<Real>(); },
          "Transmitter ECI state [r; v] [m, m/s] at transmit epoch")
      .def_readwrite("tx_clock_bias_s", &GnssChannel::tx_clock_bias_s, "Transmitter clock bias [s]")
      .def_readwrite("shapiro_delay_m", &GnssChannel::shapiro_delay_m,
                     "Shapiro (gravitational) propagation delay [m]")
      .def_readwrite("cn0_dbhz", &GnssChannel::cn0_dbhz, "Carrier-to-noise density ratio [dBHz]")
      .def_readwrite("sigma_pseudorange_m", &GnssChannel::sigma_pseudorange_m,
                     "Pseudorange measurement noise standard deviation [m]")
      .def_readwrite("sigma_doppler_hz", &GnssChannel::sigma_doppler_hz,
                     "Doppler measurement noise standard deviation [Hz]")
      .def_readwrite("sigma_carrier_phase_cycles", &GnssChannel::sigma_carrier_phase_cycles,
                     "Carrier-phase measurement noise standard deviation [cycles]");

  // ---- GNSSMeasurementsEpoch -------------------------------------------------

  py::class_<GNSSMeasurementsEpoch>(m, "GNSSMeasurementsEpoch",
                                    "GNSS measurements at a single receive epoch: the visible "
                                    "channels plus their stacked observable values/Jacobian")
      .def(py::init<>(), "Construct an empty measurement epoch")
      .def_readwrite("receive_time", &GNSSMeasurementsEpoch::receive_time,
                     "Receiver signal-reception epoch [s] in receive_time_scale")
      .def_readwrite("channels", &GNSSMeasurementsEpoch::channels,
                     "Visible/usable GNSS channels at this epoch");

  // ---- GnssConstellation -----------------------------------------------------

  py::class_<GnssConstellation, std::shared_ptr<GnssConstellation>>(
      m, "GnssConstellation",
      "Constellation of GNSS satellites with precomputed ephemerides; provides satellite "
      "state interpolation and transmitter antenna/power metadata")
      .def(py::init<>(), "Construct an empty constellation (no PRNs/ephemerides set)")
      .def(py::init<GnssConst>(), py::arg("gnss_const"),
           "Construct an empty constellation for a given GNSS system")
      .def(
          "set_satellite_states",
          [](GnssConstellation& gc, const std::vector<int>& prns, const VecXd& t_tai,
             const std::vector<MatXd>& rv_eci) {
            // rv_eci[i] is [M×6]; the C++ API expects the same layout.
            std::vector<MatXd> rv_real;
            rv_real.reserve(rv_eci.size());
            for (const auto& m : rv_eci) rv_real.push_back(m.cast<Real>());
            gc.SetSatelliteStates(prns, t_tai.cast<Real>(), rv_real);
          },
          py::arg("prns"), py::arg("t_tai"), py::arg("rv_eci"),
          "Set precomputed satellite ECI position/velocity histories (M×6 per PRN)")
      .def(
          "setup_satellite_states_from_files",
          [](GnssConstellation& gc, const std::vector<std::string>& sp3_paths,
             const std::string& antex_path, const VecXd& t_tai, GnssFreq freq,
             const std::vector<int>& prns) {
            gc.SetupSatelliteStatesFromFiles(ToPaths(sp3_paths), std::filesystem::path(antex_path),
                                             t_tai.cast<Real>(), freq, prns);
          },
          py::arg("sp3_paths"), py::arg("antex_path"), py::arg("t_tai"),
          py::arg("freq") = GnssFreq::L1, py::arg("prns") = std::vector<int>{},
          "Build ephemerides from SP3 + ANTEX files with PCO correction")
      .def(
          "load_ephemeris",
          [](GnssConstellation& gc, const std::string& filepath) {
            gc.LoadEphemeris(std::filesystem::path(filepath));
          },
          py::arg("filepath"), "Load satellite ephemerides from HDF5 file")
      .def(
          "save_ephemeris",
          [](const GnssConstellation& gc, const std::string& filepath) {
            gc.SaveEphemeris(std::filesystem::path(filepath));
          },
          py::arg("filepath"), "Save satellite ephemerides to HDF5 file")
      .def("setup_transmitters", &GnssConstellation::SetupTransmitters,
           "Load transmitter antenna patterns and power for all PRNs (GPS/Galileo/QZSS)")
      .def("get_num_satellites", &GnssConstellation::GetNumSatellites,
           "Number of satellites (PRNs) currently configured")
      .def("get_prns", &GnssConstellation::GetPrns, "List of PRNs currently configured")
      .def("get_gnss_const", &GnssConstellation::GetGnssConst,
           "GNSS system (GPS, Galileo, BDS, QZSS, ...) of this constellation")
      .def("set_fault_prns", &GnssConstellation::SetFaultPrns, py::arg("prns"),
           "Mark the given PRNs as faulted so they are reported unavailable")
      .def("is_fault_prn", &GnssConstellation::IsFaultPrn, py::arg("prn"),
           "True if `prn` has been marked faulted via set_fault_prns")
      .def(
          "get_satellite_state_eci",
          [](const GnssConstellation& gc, int prn, double t_tai) -> VecXd {
            return gc.GetSatelliteStateEci(prn, static_cast<Real>(t_tai)).cast<double>();
          },
          py::arg("prn"), py::arg("t_tai"),
          "Interpolated ECI [r; v] [m, m/s] of satellite `prn` at `t_tai` (TAI seconds)")
      .def("has_transmitter_info", &GnssConstellation::HasTransmitterInfo, py::arg("prn"),
           py::arg("freq"),
           "True if antenna and transmit-power metadata exist for `prn` and `freq`")
      .def("get_transmitter_antenna", &GnssConstellation::GetTransmitterAntenna, py::arg("prn"),
           py::arg("freq"), py::return_value_policy::reference_internal,
           "Transmit antenna pattern for `prn` and `freq`")
      .def(
          "get_transmit_power_dbw",
          [](const GnssConstellation& gc, int prn, GnssFreq freq) -> double {
            return static_cast<double>(gc.GetTransmitPowerDbw(prn, freq));
          },
          py::arg("prn"), py::arg("freq"), "Transmit power [dB-W] for `prn` and `freq`");

  // ---- GNSSMeasurements ------------------------------------------------------

  py::class_<GNSSMeasurements>(
      m, "GNSSMeasurements",
      "Receiver-side GNSS measurement generator: builds visible channels "
      "and computes pseudorange/Doppler/carrier-phase observables from one "
      "or more constellations")
      .def(py::init([](std::shared_ptr<GnssConstellation> constellation) {
             return GNSSMeasurements(constellation);
           }),
           py::arg("constellation"), "Construct with a single GNSS constellation on L1")
      .def("add_constellation", &GNSSMeasurements::AddConstellation, py::arg("constellation"),
           py::arg("frequency"),
           "Append a (constellation, frequency) pair; BuildChannels merges channels from all pairs")
      .def("set_frequency", &GNSSMeasurements::SetFrequency, py::arg("frequency"),
           "Set the carrier frequency for all current constellations (and default for future ones)")
      .def("set_options", &GNSSMeasurements::SetOptions, py::arg("options"),
           "Replace the measurement options used by build_channels()/compute()")
      .def("get_options", &GNSSMeasurements::GetOptions,
           py::return_value_policy::reference_internal, "Measurement options currently in use")
      .def("set_occluding_bodies", &GNSSMeasurements::SetOccludingBodies, py::arg("bodies"),
           "Set the spherical bodies checked for line-of-sight occlusion during visibility")
      .def("set_receiver_params", &GNSSMeasurements::SetReceiverParams, py::arg("params"),
           "Set the receiver tracking-loop/link-budget parameters used to derive noise sigmas")
      .def("set_receiver_antenna", &GNSSMeasurements::SetReceiverAntenna, py::arg("antenna"),
           "Set the receiver antenna gain pattern used in the CN0 link budget")
      .def("set_cn0_threshold", &GNSSMeasurements::SetCN0Threshold, py::arg("cn0_threshold_dbhz"),
           "Set both acquisition and tracking CN0 thresholds to the same value [dBHz]")
      .def("reset_tracking", &GNSSMeasurements::ResetTracking,
           "Clear the internal tracking state; all satellites must re-acquire on next call")
      .def(
          "set_sun_position_provider",
          [](GNSSMeasurements& meas, py::object fn) {
            meas.SetSunPositionProvider([fn](Real t) -> Vec3 {
              py::gil_scoped_acquire gil;
              return fn(static_cast<double>(t)).cast<Vec3>();
            });
          },
          py::arg("fn"),
          "Set callback returning Sun position [m] in options.frame at epoch t [TAI s]")
      .def(
          "set_boresight_target_provider",
          [](GNSSMeasurements& meas, py::object fn) {
            meas.SetBoresightTargetProvider([fn](Real t) -> Vec3 {
              py::gil_scoped_acquire gil;
              return fn(static_cast<double>(t)).cast<Vec3>();
            });
          },
          py::arg("fn"),
          "Set callback returning receiver boresight target position [m] at epoch t [TAI s]")
      .def(
          "build_channels",
          [](GNSSMeasurements& meas, double t, const VecXd& state) {
            State s(static_cast<int>(state.size()));
            s = state.cast<Real>();
            return meas.BuildChannels(static_cast<Real>(t), s);
          },
          py::arg("t"), py::arg("state"),
          "Build visible GNSS channels at receive epoch `t` [TAI s] for receiver `state`")
      .def(
          "compute",
          [](GNSSMeasurements& meas, double t, const VecXd& state) {
            State s(static_cast<int>(state.size()));
            s = state.cast<Real>();
            return meas.Compute(static_cast<Real>(t), s);
          },
          py::arg("t"), py::arg("state"),
          "Compute GNSS measurement epoch (channels + observables) at `t` for `state`")
      .def(
          "precompute",
          [](GNSSMeasurements& meas, const std::vector<double>& times,
             const std::vector<VecXd>& states, bool compute_jacobians) {
            std::vector<Real> ts(times.begin(), times.end());
            std::vector<State> ss;
            ss.reserve(states.size());
            for (const auto& v : states) {
              State s(static_cast<int>(v.size()));
              s = v.cast<Real>();
              ss.push_back(std::move(s));
            }
            return meas.Precompute(ts, ss, compute_jacobians);
          },
          py::arg("times"), py::arg("states"), py::arg("compute_jacobians") = false,
          "Batch-compute GNSS measurement epochs for a series of receive times / states");
}
