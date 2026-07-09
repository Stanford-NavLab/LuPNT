/**
 * @file antex_loader.h
 * @author Stanford NAV LAB
 * @brief ANTEX (antenna exchange format) file loader
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 * Parses IGS ANTEX (`.atx`) files and provides satellite antenna
 * phase-center-offset (PCO) queries. This is the C++ counterpart of
 * `pylupnt.interfaces.antex_file_loader.ANTEXLoader`
 * (`python/pylupnt/interfaces/antex_file_loader.py`); it focuses on the PCO
 * extraction & application workflow demonstrated in
 * `projects/Plasmasphere_Delay_Datagen/generate_antenna_pco.ipynb` (NOAZI
 * PCV pattern parsing is intentionally omitted -- only PCO is needed to set
 * up the antenna transmit positions).
 */
#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/devices/gnss_device.h"

namespace lupnt {

  /// @brief Loader / querier for satellite antenna phase-center offsets (PCO)
  /// from an IGS ANTEX file.
  ///
  /// PCO values are stored (and returned) in the satellite-antenna
  /// North/East/Up (NEU) frame, in meters (ANTEX stores them in millimeters).
  /// Converting a PCO to a geocentric/inertial frame requires the satellite
  /// "IJK" rotation matrix; see `ComputeIjkToEcefRotation` and
  /// `ApplyPcoCorrectionEcef` below, which port the convention used by
  /// `projects/GNSSEphemeris/phase_center_offset.py::ijk_to_ecef_rot`
  /// (verified against `projects/.../generate_antenna_pco.ipynb`).
  class AntexLoader {
  public:
    AntexLoader() = default;

    /// @brief Construct and load a single ANTEX file.
    explicit AntexLoader(const std::filesystem::path& filepath);

    /// @brief Parse an additional ANTEX file and merge its satellite entries in.
    void LoadFile(const std::filesystem::path& filepath);

    /// @brief Look up the satellite antenna phase-center offset (PCO), in the
    /// satellite North/East/Up (NEU) antenna frame [m], for satellite `prn`
    /// of constellation `gnss_const`, frequency `freq`, at epoch `t_tai`
    /// (used to select the correct validity-period entry; ANTEX validity
    /// spans are wide -- weeks to years -- so any reasonable choice of time
    /// system for `t_tai` is adequate). Mirrors `ANTEXLoader.get_pco`.
    /// Throws if no matching entry/frequency is found.
    Vec3d GetPco(GnssConst gnss_const, int prn, GnssFreq freq, Real t_tai) const;

    /// @brief Same as `GetPco`, but specifying the raw ANTEX frequency code
    /// directly (e.g. `"G01"`, `"E05"`) instead of a `GnssFreq` enum value.
    Vec3d GetPco(const std::string& gnss_letter, int prn, const std::string& antex_freq_code,
                 Real t_tai) const;

    /// @brief True if a PCO is available for this satellite/frequency at `t_tai`.
    /// ANTEX files legitimately omit frequencies a satellite does not transmit (e.g. L5 on
    /// older GPS blocks), so callers can fall back to a zero offset instead of throwing.
    bool HasPco(GnssConst gnss_const, int prn, GnssFreq freq, Real t_tai) const;

    /// @brief List the ANTEX frequency codes available for satellite `prn`
    /// of constellation `gnss_const` at epoch `t_tai` (e.g. `{"G01","G02","G05"}`).
    std::vector<std::string> GetAvailableFreqCodes(GnssConst gnss_const, int prn, Real t_tai) const;

    bool HasSatellite(GnssConst gnss_const, int prn) const;

    // ---- PCO -> ECEF correction --------------------------------------------

    /// @brief "IJK" rotation matrix used to express a satellite-antenna NEU
    /// PCO vector in ECEF, as `pco_ecef = Cijk * pco_neu`. Columns are:
    ///   ivec = jvec x kvec, jvec = normalize(r_sun_ecef - r_sat_ecef),
    ///   kvec = -normalize(r_sat_ecef)
    /// **Note**: this triad is intentionally distinct from (and, in general,
    /// not orthonormal with) the canonical `GnssAttitude` body frame --
    /// `jvec` here is the *raw* (un-orthogonalized) approximate Earth-Sun
    /// direction, not `GnssAttitude::GetEy()`. It is ported verbatim from
    /// `phase_center_offset.py::ijk_to_ecef_rot` (the reference used to
    /// generate the precomputed PCO-corrected ephemerides in
    /// `generate_antenna_pco.ipynb` / `test_ephemeris_loading.ipynb`), and
    /// MUST be kept distinct to remain numerically consistent with that data.
    static Mat3d ComputeIjkToEcefRotation(Real t_tai, const Vec3d& r_sat_ecef);

    /// @brief Apply the PCO correction to an SP3 (center-of-mass) ECEF
    /// position: `pos_corrected = pos_sp3_ecef + Cijk(t_tai, pos_sp3_ecef) * pco_neu`.
    /// Mirrors the workflow in `generate_antenna_pco.ipynb`.
    static Vec3d ApplyPcoCorrectionEcef(Real t_tai, const Vec3d& pos_sp3_ecef,
                                        const Vec3d& pco_neu_m);

    // ---- Identifier helpers (also used by `GnssConstellation` to build SP3 /
    // ANTEX satellite identifiers, e.g. "G05", from `(GnssConst, prn)`) -------

    /// @brief Single-letter GNSS system identifier used by SP3/ANTEX/RINEX
    /// satellite IDs, e.g. `GnssConst::GPS` -> `"G"`, `GnssConst::GALILEO` -> `"E"`.
    static std::string GnssLetter(GnssConst gnss_const);

    /// @brief SP3/ANTEX/RINEX satellite identifier, e.g. `(GPS, 5)` -> `"G05"`.
    static std::string SatId(GnssConst gnss_const, int prn);

  private:
    struct FreqPattern {
      std::string freq_code;
      bool has_pco = false;
      Vec3d pco_neu_m = Vec3d::Zero();  // [N, E, U] meters
    };

    struct SatAntennaEntry {
      std::string sat_id;  // e.g. "G05"
      bool has_valid_from = false;
      bool has_valid_until = false;
      double valid_from_tai = 0.0;
      double valid_until_tai = 0.0;
      std::map<std::string, FreqPattern> freqs;  // ANTEX freq code -> pattern
    };

    std::map<std::string, std::vector<SatAntennaEntry>> sat_index_;  // sat_id -> entries

    void ParseFile(const std::filesystem::path& filepath);

    static std::string FreqToAntexCode(const std::string& gnss_letter,
                                       const std::string& freq_name);

    const SatAntennaEntry& SelectEntry(const std::string& sat_id, double t_tai) const;
  };

}  // namespace lupnt
