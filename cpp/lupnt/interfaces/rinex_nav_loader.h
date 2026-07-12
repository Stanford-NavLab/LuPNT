/**
 * @file rinex_nav_loader.h
 * @author Stanford NAV LAB
 * @brief RINEX navigation (broadcast ephemeris / BRDC) file loader
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 * Parses RINEX V3 navigation ("BRDC") files and provides broadcast-ephemeris
 * ECEF position/velocity/clock-correction queries via Keplerian propagation
 * of the navigation message nearest the requested epoch. This is the C++
 * counterpart of `pylupnt.interfaces.gnss_file_loader.BRDCLoader`
 * (`python/pylupnt/interfaces/gnss_file_loader.py`); it supports the
 * Keplerian-element systems GPS / Galileo / BeiDou / QZSS (`G`/`E`/`C`/`J`)
 * -- GLONASS (`R`, which uses tabulated orbital state vectors rather than
 * Keplerian elements) is intentionally not supported, matching
 * `BRDCLoader.get_posvelclock_all`'s default exclusion of GLONASS.
 *
 * Simplification vs. the Python reference: the Galileo-specific GST/GPST
 * "GAGP" system-time-correction term (`t_corr_sys`, parsed from the
 * "TIME SYSTEM CORR" header record) is omitted (`t_corr_sys = 0`). This term
 * only affects the returned satellite *clock* correction (sub-100ns, i.e.
 * sub-30m range-equivalent) and has **no effect whatsoever on the broadcast
 * position/velocity** -- which is this loader's primary purpose (setting up
 * antenna transmit positions / link-budget geometry).
 */
#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Loader / propagator for RINEX V3 broadcast-ephemeris ("BRDC")
  /// navigation files.
  ///
  /// Satellites are identified the same way as in `Sp3Loader`, e.g. `"G01"`.
  class RinexNavLoader {
  public:
    RinexNavLoader() = default;

    /// @brief Construct and load a single RINEX nav file.
    explicit RinexNavLoader(const std::filesystem::path& filepath);

    /// @brief Construct and load multiple RINEX nav files (e.g. consecutive
    /// days); navigation messages are concatenated per-satellite.
    explicit RinexNavLoader(const std::vector<std::filesystem::path>& filepaths);

    /// @brief Parse an additional RINEX nav file and merge its messages in.
    void LoadFile(const std::filesystem::path& filepath);

    /// @brief Parse a YUMA-format GPS almanac (e.g. from CelesTrak / USCG NAVCEN, a small
    /// no-authentication text download) and merge its per-PRN orbital slots in as broadcast
    /// navigation messages. YUMA carries only the coarse Keplerian set (e, toa, i, OMEGADOT,
    /// sqrt_a, OMEGA0, omega, M0, af0, af1, week); the higher-order harmonic corrections and
    /// delta_n / idot are set to zero. This is intended as an orbital-slot *seed* for numerical
    /// propagation (see the lunar-GNSS ODTS `almanac` constellation source), not as a precise
    /// broadcast product. Only GPS (`G`) satellites are produced.
    void LoadYumaFile(const std::filesystem::path& filepath);

    /// @brief Reference epoch [TAI seconds] of the most recent navigation message loaded for
    /// `sat_id` (the freshest orbital slot). Used to pick the seed epoch at which the broadcast
    /// Keplerian elements are converted to a Cartesian state for numerical propagation, so the
    /// state is taken at the message's own time-of-ephemeris (t_k = 0) rather than
    /// Keplerian-extrapolated far from it.
    double GetLatestEpochTai(const std::string& sat_id) const;

    /// @brief Build the CDDIS BRDC filename covering `epoch`, e.g.
    /// `BRDC00IGS_R_20260140000_01D_MN.rnx` (uncompressed).
    static std::string FilenameForEpoch(Real epoch, Time time_scale = Time::UTC);

    /// @brief Download/cache the multi-GNSS broadcast-ephemeris ("BRDC") RINEX
    /// nav file covering `epoch`, then return the uncompressed `.rnx` path.
    ///
    /// The default cache directory is `GetOutputDir("gnss_files") / "brdc"`.
    /// Mirrors `Sp3Loader::DownloadFileForEpoch`: downloads use
    /// `curl -L --netrc-optional`, so credentials can be supplied via `~/.netrc`
    /// or the `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` environment variables.
    /// The daily product is tried at both CDDIS layouts
    /// (`.../daily/YYYY/DDD/YYp/` then `.../daily/YYYY/brdc/`). Throws if CDDIS
    /// returns an Earthdata Login HTML page.
    static std::filesystem::path DownloadFileForEpoch(Real epoch, Time time_scale = Time::UTC,
                                                      const std::filesystem::path& cache_dir
                                                      = std::filesystem::path());

    /// @brief Identifiers of all satellites with navigation messages loaded
    /// so far (e.g. `{"G01", "G02", ..., "E11", ...}`); GLONASS excluded.
    const std::vector<std::string>& GetSatellites() const { return sats_; }

    bool HasSatellite(const std::string& sat_id) const;

    /// @brief Broadcast-ephemeris ECEF position/velocity [m, m/s] (`rv_ecef`)
    /// and clock correction [s] (`clock_corr_s`, polynomial + relativistic
    /// terms only -- see file-level note on the omitted Galileo system-time
    /// term) of satellite `sat_id` at `t_tai`, computed via Keplerian
    /// propagation of the navigation message with the closest time-of-epoch.
    /// Mirrors `BRDCLoader.get_posvelclock` for systems G/E/C/J.
    void GetPosVelClock(const std::string& sat_id, Real t_tai, Vec6& rv_ecef,
                        Real& clock_corr_s) const;

    /// @brief Broadcast-ephemeris ECEF position/velocity [m, m/s] only.
    Vec6 GetPosVel(const std::string& sat_id, Real t_tai) const;

  private:
    /// @brief A single broadcast navigation message (Keplerian orbital
    /// elements + clock polynomial), for systems G/E/C/J. Field names follow
    /// RINEX 3.04 (see http://acc.igs.org/misc/rinex304.pdf).
    struct NavMessage {
      double epoch_tai = 0.0;
      double af0 = 0.0, af1 = 0.0, af2 = 0.0;
      double crs = 0.0, delta_n = 0.0, m0 = 0.0;
      double cuc = 0.0, ecc = 0.0, cus = 0.0, sqrt_a = 0.0;
      double toe = 0.0, cic = 0.0, omega0 = 0.0, cis = 0.0;
      double i0 = 0.0, crc = 0.0, omega = 0.0, omega_dot = 0.0, idot = 0.0;
      double week = 0.0;
    };

    std::vector<std::string> sats_;
    std::map<std::string, std::vector<NavMessage>> nav_;  // sat_id -> messages (any order)

    void ParseFile(const std::filesystem::path& filepath);

    /// @brief Index of the navigation message in `nav_[sat_id]` whose
    /// `epoch_tai` is closest to `t_tai` (mirrors `np.argmin(abs(epoch_tai - t))`).
    int FindClosestMessage(const std::string& sat_id, double t_tai) const;
  };

}  // namespace lupnt
