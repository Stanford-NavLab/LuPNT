/**
 * @file sp3_loader.h
 * @author Stanford NAV LAB
 * @brief SP3 (precise GNSS ephemeris) file loader
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 * Parses IGS SP3 precise-ephemeris files and provides interpolated
 * ECEF position / velocity / clock-bias queries. This is the C++ counterpart
 * of `pylupnt.interfaces.gnss_file_loader.SP3Loader`
 * (`python/pylupnt/interfaces/gnss_file_loader.py`); unlike the Python
 * version can download/cache modern CDDIS SP3 files for requested epochs, or
 * parse files that are already present on disk.
 *
 * Position/velocity interpolation is provided by a per-satellite piecewise
 * Chebyshev fit over the tabulated ECEF position and clock samples. Query-time
 * velocity is the analytic derivative of the fitted position polynomial.
 */
#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/numerics/cheby_fit.h"

namespace lupnt {

  /// @brief Loader / interpolator for IGS SP3 precise-ephemeris files.
  ///
  /// Satellites are identified by their SP3 identifiers, e.g. `"G01"` (GPS
  /// PRN 1), `"E11"` (Galileo PRN 11), `"C06"` (BeiDou PRN 6), `"J03"`
  /// (QZSS PRN 3).
  class Sp3Loader {
  public:
    Sp3Loader() = default;

    /// @brief Construct and load a single SP3 file.
    explicit Sp3Loader(const std::filesystem::path& filepath);

    /// @brief Construct and load multiple SP3 files (e.g. consecutive days);
    /// samples are concatenated and sorted by epoch per satellite.
    explicit Sp3Loader(const std::vector<std::filesystem::path>& filepaths);

    /// @brief Parse an additional SP3 file and merge its samples in.
    void LoadFile(const std::filesystem::path& filepath);

    /// @brief Build the modern CDDIS/COD MGEX final SP3 filename for `epoch`.
    ///
    /// @param epoch      Epoch value in `time_scale` seconds.
    /// @param time_scale Time scale of `epoch` (UTC by default).
    /// @return Filename such as `COD0MGXFIN_20250010000_01D_05M_ORB.SP3`.
    static std::string FilenameForEpoch(Real epoch, Time time_scale = Time::UTC);

    /// @brief Build the CDDIS download URL for the SP3 product covering `epoch`.
    ///
    /// This follows the same modern COD/MGEX final-product convention as the
    /// Python `SP3Loader`. Pre-2017 legacy `.Z` products are intentionally not
    /// supported by the C++ helper.
    static std::string UrlForEpoch(Real epoch, Time time_scale = Time::UTC);

    /// @brief Download/cache the SP3 file covering `epoch`, then return the
    /// uncompressed `.SP3` path.
    ///
    /// The default cache directory is `GetOutputDir("gnss_files") / "sp3"`,
    /// matching the Python loader. Downloads use `curl -L --netrc-optional`,
    /// so credentials can be supplied with `~/.netrc`; if
    /// `EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD` are set, they are also
    /// passed to curl. Throws if CDDIS returns an Earthdata Login HTML page.
    static std::filesystem::path DownloadFileForEpoch(Real epoch, Time time_scale = Time::UTC,
                                                      const std::filesystem::path& cache_dir
                                                      = std::filesystem::path());

    /// @brief SP3 identifiers of all satellites with data loaded so far,
    /// e.g. `{"G01", "G02", ..., "E11", ...}`.
    const std::vector<std::string>& GetSatellites() const { return sats_; }

    bool HasSatellite(const std::string& sat_id) const;

    /// @brief Interpolated ECEF position [m] / velocity [m/s] (`rv_ecef`,
    /// columns 0-2 / 3-5) and clock bias [s] (`clock_bias_s`, including the
    /// relativistic correction `-2/c^2 * dot(r, v)`, mirroring
    /// `SP3Loader.get_posvelclock`) of satellite `sat_id` at epoch `t_tai`
    /// (TAI seconds). Throws if `sat_id` is unknown or `t_tai` is outside
    /// the loaded ephemeris span (no orbital-propagation fallback, unlike
    /// the Python reference -- the precise-ephemeris window is expected to
    /// cover the requested epochs).
    void GetPosVelClock(const std::string& sat_id, Real t_tai, Vec6& rv_ecef,
                        Real& clock_bias_s) const;

    /// @brief Interpolated ECEF position [m] / velocity [m/s] only.
    Vec6 GetPosVel(const std::string& sat_id, Real t_tai) const;

    /// @brief Time span [TAI seconds] covered by the loaded ephemeris for
    /// satellite `sat_id` (`{t_min, t_max}`).
    std::pair<double, double> GetTimeSpan(const std::string& sat_id) const;

  private:
    std::vector<std::string> sats_;            // SP3 identifiers, e.g. "G01"
    std::map<std::string, VecXd> epochs_tai_;  // sat -> [N] epochs (TAI seconds, sorted)
    std::map<std::string, MatXd> pos_clock_;   // sat -> [N x 4] ECEF (x,y,z) [m], clock [s]
    std::map<std::string, ChebyshevFitModel> pos_clock_cheby_;  // sat -> [x,y,z,clock] fit

    void ParseFile(const std::filesystem::path& filepath);
    void RebuildChebyshevModels();
  };

}  // namespace lupnt
