
#include "lupnt/data/tai_utc.h"

#include <lupnt/core/constants.h>
#include <lupnt/core/file.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace lupnt {

  std::mutex tai_utc_mutex;
  Ptr<TaiUtcFileData> tai_utc_data;

  void LoadTaiUtcFileData(const std::filesystem::path& filepath) {
    std::lock_guard<std::mutex> lock(tai_utc_mutex);
    if (tai_utc_data) return;  // Data already loaded

    size_t n_lines = CountLines(filepath.string());
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    // Initialize TaiUtcFileData struct
    tai_utc_data = MakePtr<TaiUtcFileData>();
    tai_utc_data->jd.resize(n_lines);
    tai_utc_data->tai_utc.resize(n_lines);
    tai_utc_data->mjd0.resize(n_lines);
    tai_utc_data->scale.resize(n_lines);

    std::string line;
    size_t index = 0;
    while (std::getline(file, line)) {
      if (line.empty()) continue;

      double jd, tai_utc, mjd0, scale;
      std::istringstream iss(line);
      std::string year, month, day, eq_jd, tai_utc_eq, s1, plus, paren_mjd, minus, paren, times, s2;

      iss >> year >> month >> day >> eq_jd >> jd >> tai_utc_eq >> tai_utc >> s1 >> plus >> paren_mjd
          >> minus >> mjd0 >> paren >> times >> scale >> s2;

      tai_utc_data->jd(index) = jd;
      tai_utc_data->tai_utc(index) = tai_utc;
      tai_utc_data->mjd0(index) = mjd0;
      tai_utc_data->scale(index) = scale;
      ++index;
    }

    file.close();
    return;
  }

  double Evaluate(double mjd, int i) {
    return tai_utc_data->tai_utc(i) + tai_utc_data->scale(i) * (mjd - tai_utc_data->mjd0(i));
  }

  double GetTaiUtcDifference(double mjd) {
    if (!tai_utc_data) LoadTaiUtcFileData(GetFilePath(TAI_UTC_FILENAME));
    double jd = mjd + JD_MJD_OFFSET;
    if (jd < tai_utc_data->jd(0)) {
      return 0.0;
    } else if (jd > tai_utc_data->jd(tai_utc_data->jd.size() - 1)) {
      int index = tai_utc_data->jd.size() - 1;
      return Evaluate(mjd, index);
    }

    // Find closest MJD for the immediate smaller date
    double* start = tai_utc_data->jd.data();
    double* end = tai_utc_data->jd.data() + tai_utc_data->jd.size();
    auto it = std::upper_bound(start, end, jd);

    int index = (it - tai_utc_data->jd.data()) - 1;
    if (index < 0) index = 0;
    return Evaluate(mjd, index);
  }
}  // namespace lupnt
