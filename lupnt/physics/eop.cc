#include "eop.h"

#include <lupnt/core/constants.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

namespace lupnt {

// Global shared pointer and mutex
std::shared_ptr<EOPData> eop_data;
std::mutex eop_mutex;

size_t CountLines(const std::string& filename) {
  std::ifstream file(filename);
  assert(file.is_open() && "Unable to open file");

  size_t lineCount = 0;
  std::string line;
  while (std::getline(file, line)) {
    ++lineCount;
  }
  file.close();
  return lineCount;
}

std::shared_ptr<EOPData> LoadEOPData(const std::string& filename) {
  std::lock_guard<std::mutex> lock(eop_mutex);
  if (eop_data) {
    return eop_data;  // Return already loaded data
  }
  int n_header_lines = 14;
  size_t n_lines = CountLines(filename) - n_header_lines;
  std::ifstream file(filename);
  assert(file.is_open() && "Unable to open file");

  // Skip header lines
  std::string line;
  for (int i = 0; i < n_header_lines; ++i) {
    std::getline(file, line);
  }

  // Initialize EOPData struct
  auto data = std::make_shared<EOPData>();
  data->years.resize(n_lines);
  data->months.resize(n_lines);
  data->days.resize(n_lines);
  data->mjds.resize(n_lines);
  data->x.resize(n_lines);
  data->y.resize(n_lines);
  data->ut1_utc.resize(n_lines);
  data->lod.resize(n_lines);
  data->dPsi.resize(n_lines);
  data->dEps.resize(n_lines);
  data->xErr.resize(n_lines);
  data->yErr.resize(n_lines);
  data->ut1_utcErr.resize(n_lines);
  data->lodErr.resize(n_lines);
  data->dPsiErr.resize(n_lines);
  data->dEpsErr.resize(n_lines);

  size_t row = 0;
  // Read data lines
  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    int year, month, day, mjd;
    double x, y, ut1_utc, lod, dPsi, dEps, xErr, yErr, ut1_utcErr, lodErr,
        dPsiErr, dEpsErr;

    iss >> year >> month >> day >> mjd >> x >> y >> ut1_utc >> lod >> dPsi >>
        dEps >> xErr >> yErr >> ut1_utcErr >> lodErr >> dPsiErr >> dEpsErr;

    data->years(row) = year;
    data->months(row) = month;
    data->days(row) = day;
    data->mjds(row) = mjd;
    data->x(row) = x;
    data->y(row) = y;
    data->ut1_utc(row) = ut1_utc;
    data->lod(row) = lod;
    data->dPsi(row) = dPsi;
    data->dEps(row) = dEps;
    data->xErr(row) = xErr;
    data->yErr(row) = yErr;
    data->ut1_utcErr(row) = ut1_utcErr;
    data->lodErr(row) = lodErr;
    data->dPsiErr(row) = dPsiErr;
    data->dEpsErr(row) = dEpsErr;

    ++row;
  }

  file.close();
  eop_data = data;
  return eop_data;
}

EOPResult InterpolateEOPData(const std::shared_ptr<EOPData>& eop_data,
                             Real mjdUTC, bool interpolate) {
  if (!eop_data) {
    std::shared_ptr<EOPData> eop_data =
        LoadEOPData(GetFilePath("eopc04_08.62-now"));
  }
  EOPResult result;

  assert(mjdUTC >= eop_data->mjds(0) &&
         mjdUTC <= eop_data->mjds(eop_data->mjds.size() - 1) &&
         "Requested MJD is outside the range of the EOP data");

  // Find closest MJD
  int* start = eop_data->mjds.data();
  int* end = eop_data->mjds.data() + eop_data->mjds.size();
  int val = static_cast<int>(mjdUTC);
  auto it = std::lower_bound(start, end, val);

  int i_prev = it - eop_data->mjds.data();
  int i_next = i_prev + 1;

  if (interpolate) {
    // Linear interpolation
    Real s = (mjdUTC - eop_data->mjds(i_prev)) /
             (eop_data->mjds(i_next) - eop_data->mjds(i_prev));

    auto interp = [](Real x0, Real x1, Real s) { return x0 + (x1 - x0) * s; };
    result.x_pole =
        interp(eop_data->x(i_prev), eop_data->x(i_next), s) * RAD_ARCSEC;
    result.y_pole =
        interp(eop_data->y(i_prev), eop_data->y(i_next), s) * RAD_ARCSEC;
    result.UT1_UTC =
        interp(eop_data->ut1_utc(i_prev), eop_data->ut1_utc(i_next), s);
    result.LOD = interp(eop_data->lod(i_prev), eop_data->lod(i_next), s);
    result.dPsi =
        interp(eop_data->dPsi(i_prev), eop_data->dPsi(i_next), s) * RAD_ARCSEC;
    result.dEps =
        interp(eop_data->dEps(i_prev), eop_data->dEps(i_next), s) * RAD_ARCSEC;
    result.dx_pole =
        interp(eop_data->xErr(i_prev), eop_data->xErr(i_next), s) * RAD_ARCSEC;
    result.dy_pole =
        interp(eop_data->yErr(i_prev), eop_data->yErr(i_next), s) * RAD_ARCSEC;
    result.TAI_UTC = eop_data->ut1_utcErr(i_prev);

  } else {
    // Use the nearest point
    int index = abs(mjdUTC - eop_data->mjds(i_prev)) <
                        abs(mjdUTC - eop_data->mjds(i_next))
                    ? i_prev
                    : i_next;
    result.x_pole = eop_data->x(index) * RAD_ARCSEC;
    result.y_pole = eop_data->y(index) * RAD_ARCSEC;
    result.UT1_UTC = eop_data->ut1_utc(index);
    result.LOD = eop_data->lod(index);
    result.dPsi = eop_data->dPsi(index) * RAD_ARCSEC;
    result.dEps = eop_data->dEps(index) * RAD_ARCSEC;
    result.dx_pole = eop_data->xErr(index) * RAD_ARCSEC;
    result.dy_pole = eop_data->yErr(index) * RAD_ARCSEC;
    result.TAI_UTC = eop_data->ut1_utcErr(index);
  }

  return result;
}

}  // namespace lupnt