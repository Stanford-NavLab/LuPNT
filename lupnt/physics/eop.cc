#include "eop.h"

#include <lupnt/core/constants.h>
#include <lupnt/core/file.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

namespace lupnt {

std::shared_ptr<EopFileData> eop;
std::mutex eop_mutex;

void LoadEopData(const std::filesystem::path& filepath) {
  std::lock_guard<std::mutex> lock(eop_mutex);
  if (eop) return;  // Data already loaded

  int n_header_lines = 14;
  size_t n_lines = CountLines(filepath) - n_header_lines;
  std::ifstream file(filepath);
  assert(file.is_open() && "Unable to open file");

  // Skip header lines
  std::string line;
  for (int i = 0; i < n_header_lines; ++i) {
    std::getline(file, line);
  }

  // Initialize EopFileData struct
  eop = std::make_shared<EopFileData>();
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
    double x, y, ut1_utc, lod, dpsi, deps, xErr, yErr, ut1_utc_err, lod_err,
        dpsi_err, deps_err;

    iss >> year >> month >> day >> mjd >> x >> y >> ut1_utc >> lod >> dpsi >>
        deps >> xErr >> yErr >> ut1_utc_err >> lod_err >> dpsi_err >> deps_err;

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

Real GetUt1UtcDifference(Real mjd_utc, bool interpolate) {
  EopData eop = GetEopData(mjd_utc, interpolate);
  return eop.ut1_utc;
}

EopData GetEopData(Real mjd_utc, bool interpolate) {
  if (!eop) LoadEopData(GetFilePath(EOP_FILENAME));

  EopData data;

  // Check if the requested MJD is outside the range
  if (mjd_utc <= eop->mjds_utc(0) ||
      mjd_utc >= eop->mjds_utc(eop->mjds_utc.size() - 1)) {
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

  // Find closest MJD
  int* start = eop->mjds_utc.data();
  int* end = eop->mjds_utc.data() + eop->mjds_utc.size();
  int val = static_cast<int>(mjd_utc);
  auto it = std::lower_bound(start, end, val);

  int i_prev = it - eop->mjds_utc.data();
  int i_next = i_prev + 1;

  if (interpolate) {
    // Linear interpolation
    Real s = (mjd_utc - eop->mjds_utc(i_prev)) /
             (eop->mjds_utc(i_next) - eop->mjds_utc(i_prev));

    auto interp = [](Real x0, Real x1, Real s) { return x0 + (x1 - x0) * s; };
    data.x_pole = interp(eop->x(i_prev), eop->x(i_next), s) * RAD_ARCSEC;
    data.y_pole = interp(eop->y(i_prev), eop->y(i_next), s) * RAD_ARCSEC;
    data.ut1_utc = interp(eop->ut1_utc(i_prev), eop->ut1_utc(i_next), s);
    data.lod = interp(eop->lod(i_prev), eop->lod(i_next), s);
    data.dpsi = interp(eop->dpsi(i_prev), eop->dpsi(i_next), s) * RAD_ARCSEC;
    data.deps = interp(eop->deps(i_prev), eop->deps(i_next), s) * RAD_ARCSEC;
    data.dx_pole = interp(eop->xErr(i_prev), eop->xErr(i_next), s) * RAD_ARCSEC;
    data.dy_pole = interp(eop->yErr(i_prev), eop->yErr(i_next), s) * RAD_ARCSEC;
    data.tai_utc = eop->ut1_utc_err(i_prev);

  } else {
    // Use the nearest point
    int index = abs(mjd_utc - eop->mjds_utc(i_prev)) <
                        abs(mjd_utc - eop->mjds_utc(i_next))
                    ? i_prev
                    : i_next;
    data.x_pole = eop->x(index) * RAD_ARCSEC;
    data.y_pole = eop->y(index) * RAD_ARCSEC;
    data.ut1_utc = eop->ut1_utc(index);
    data.lod = eop->lod(index);
    data.dpsi = eop->dpsi(index) * RAD_ARCSEC;
    data.deps = eop->deps(index) * RAD_ARCSEC;
    data.dx_pole = eop->xErr(index) * RAD_ARCSEC;
    data.dy_pole = eop->yErr(index) * RAD_ARCSEC;
    data.tai_utc = eop->ut1_utc_err(index);
  }

  return data;
}

}  // namespace lupnt