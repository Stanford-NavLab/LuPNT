#include "lupnt/data/eop.h"

#include <lupnt/core/constants.h>
#include <lupnt/core/file.h>
#include <lupnt/numerics/interpolation.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace lupnt {

  Ptr<EopFileData> eop;
  std::mutex eop_mutex;

  void LoadEopFileData(const std::filesystem::path& filepath) {
    std::lock_guard<std::mutex> lock(eop_mutex);
    if (eop) return;  // Data already loaded

    int n_header_lines = 14;
    size_t n_lines = CountLines(filepath.string()) - n_header_lines;
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    // Skip header lines
    std::string line;
    for (int i = 0; i < n_header_lines; ++i) {
      std::getline(file, line);
    }

    // Initialize EopFileData struct
    eop = MakePtr<EopFileData>();
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

  Real GetUt1UtcDifference(Real mjd_utc) {
    EopData eop = GetEopData(mjd_utc);
    return eop.ut1_utc;
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
