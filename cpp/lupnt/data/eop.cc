#include "lupnt/data/eop.h"

#include <fmt/format.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "lupnt/core/constants.h"
#include "lupnt/core/file.h"
#include "lupnt/core/logger.h"
#include "lupnt/numerics/interpolation.h"

namespace lupnt {

  UniquePtr<EopFileData> eop;
  std::mutex eop_mutex;

  namespace {
    /// Download `url` to `dest_path` using curl. The file is first written to a
    /// temporary path and atomically moved into place on success, leaving any
    /// existing file at `dest_path` untouched on failure.
    bool DownloadFile(const std::string& url, const std::filesystem::path& dest_path) {
      std::filesystem::path tmp_path = dest_path;
      tmp_path += ".part";
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);

      std::filesystem::create_directories(dest_path.parent_path(), ec);

      std::string cmd
          = fmt::format("curl -fsSL --connect-timeout 10 --max-time 600 -o \"{}\" \"{}\"",
                        tmp_path.string(), url);
      bool ok = std::system(cmd.c_str()) == 0 && std::filesystem::exists(tmp_path)
                && std::filesystem::file_size(tmp_path) > 0;
      if (!ok) {
        std::filesystem::remove(tmp_path, ec);
        return false;
      }
      std::filesystem::rename(tmp_path, dest_path, ec);
      if (ec) {
        ec.clear();
        std::filesystem::copy_file(tmp_path, dest_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(tmp_path, ec);
      }
      return !ec;
    }
  }  // namespace

  void LoadEopFileData(const std::filesystem::path& filepath, bool force) {
    std::lock_guard<std::mutex> lock(eop_mutex);
    if (eop && !force) return;  // Data already loaded

    int n_header_lines = 14;
    size_t n_lines = CountLines(filepath.string()) - n_header_lines;
    std::ifstream file = OpenFile<std::ifstream>(filepath);
    // Skip header lines
    std::string line;
    for (int i = 0; i < n_header_lines; ++i) {
      std::getline(file, line);
    }

    // Initialize EopFileData struct
    eop = MakeUnique<EopFileData>();
    eop->years.resize(n_lines);
    eop->months.resize(n_lines);
    eop->days.resize(n_lines);
    eop->mjds_utc.resize(n_lines);
    eop->x.resize(n_lines);
    eop->y.resize(n_lines);
    eop->ut1_utc.resize(n_lines);
    eop->lod.resize(n_lines);
    eop->dpsi.resize(n_lines);
    eop->deps.resize(n_lines);
    eop->xErr.resize(n_lines);
    eop->yErr.resize(n_lines);
    eop->ut1_utc_err.resize(n_lines);
    eop->lod_err.resize(n_lines);
    eop->dpsi_err.resize(n_lines);
    eop->deps_err.resize(n_lines);

    size_t row = 0;
    // Read data lines
    while (std::getline(file, line)) {
      if (line.empty()) {
        continue;
      }
      std::istringstream iss(line);
      int year, month, day, mjd;
      double x, y, ut1_utc, lod, dpsi, deps, xErr, yErr, ut1_utc_err, lod_err, dpsi_err, deps_err;

      iss >> year >> month >> day >> mjd >> x >> y >> ut1_utc >> lod >> dpsi >> deps >> xErr >> yErr
          >> ut1_utc_err >> lod_err >> dpsi_err >> deps_err;

      eop->years(row) = year;
      eop->months(row) = month;
      eop->days(row) = day;
      eop->mjds_utc(row) = mjd;
      eop->x(row) = x;
      eop->y(row) = y;
      eop->ut1_utc(row) = ut1_utc;
      eop->lod(row) = lod;
      eop->dpsi(row) = dpsi;
      eop->deps(row) = deps;
      eop->xErr(row) = xErr;
      eop->yErr(row) = yErr;
      eop->ut1_utc_err(row) = ut1_utc_err;
      eop->lod_err(row) = lod_err;
      eop->dpsi_err(row) = dpsi_err;
      eop->deps_err(row) = deps_err;

      ++row;
    }

    file.close();
    return;
  }

  bool LoadLatestEopFromIers(bool force) {
    static const std::string url
        = "https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt";
    std::filesystem::path dest_path = GetDataPath() / "planetary_coeff" / "eopc04_iers_latest.txt";

    if (DownloadFile(url, dest_path)) {
      Logger::Info(fmt::format("EOP: downloaded latest IERS EOP data to '{}'", dest_path.string()),
                   "Eop");
      LoadEopFileData(dest_path, force);
      return true;
    }

    Logger::Warn("EOP: failed to download latest IERS EOP data, falling back to bundled data",
                 "Eop");
    GetEopFileData();  // Ensure the bundled file is loaded as a fallback
    return false;
  }

  Real GetUt1UtcDifference(Real mjd_utc) {
    EopData eop = GetEopData(mjd_utc);
    return eop.ut1_utc;
  }

  EopFileData* GetEopFileData() {
    if (!eop) LoadEopFileData(GetFilePath(EOP_FILENAME));
    return eop.get();
  }

  EopData GetEopData(Real mjd_utc) {
    if (!eop) LoadEopFileData(GetFilePath(EOP_FILENAME));

    EopData data;

    // Check if the requested MJD is outside the range
    if (mjd_utc <= eop->mjds_utc(0) || mjd_utc >= eop->mjds_utc(eop->mjds_utc.size() - 1)) {
      int i = mjd_utc <= eop->mjds_utc(0) ? 0 : eop->mjds_utc.size() - 1;
      data.x_pole = eop->x(i) * RAD_ARCSEC;
      data.y_pole = eop->y(i) * RAD_ARCSEC;
      data.ut1_utc = eop->ut1_utc(i);
      data.lod = eop->lod(i);
      data.dpsi = eop->dpsi(i) * RAD_ARCSEC;
      data.deps = eop->deps(i) * RAD_ARCSEC;
      data.dx_pole = eop->xErr(i) * RAD_ARCSEC;
      data.dy_pole = eop->yErr(i) * RAD_ARCSEC;
      data.tai_utc = eop->ut1_utc_err(i);
      return data;
    }

    int order = 3;
    LagrangeInterpolator interp(eop->mjds_utc, mjd_utc.val(), order);
    data.x_pole = interp.Interpolate(eop->x) * RAD_ARCSEC;
    data.y_pole = interp.Interpolate(eop->y) * RAD_ARCSEC;
    data.ut1_utc = interp.Interpolate(eop->ut1_utc);
    data.lod = interp.Interpolate(eop->lod);
    data.dpsi = interp.Interpolate(eop->dpsi) * RAD_ARCSEC;
    data.deps = interp.Interpolate(eop->deps) * RAD_ARCSEC;
    data.dx_pole = interp.Interpolate(eop->xErr) * RAD_ARCSEC;
    data.dy_pole = interp.Interpolate(eop->yErr) * RAD_ARCSEC;
    data.tai_utc = interp.Interpolate(eop->ut1_utc_err);
    return data;
  }

}  // namespace lupnt
