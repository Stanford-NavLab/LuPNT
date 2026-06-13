/**
 * @file string_file_utils.cpp
 * @author Keidai Iiyama
 * @brief This file contains utility functions for string and file manipulation.
 * @version 0.1
 * @date 2025-02-17
 */

#include "lupnt/environment/plasma/core/string_file_utils.h"

#include <cctype>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "lupnt/environment/plasma/core/user_filepath.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {

  std::string trim_right_copy(const std::string& input) {
    auto it = input.rbegin();
    while (it != input.rend() && std::isspace(*it)) {
      ++it;
    }
    return std::string(input.begin(), it.base());
  }

  std::vector<std::string> split_string(const std::string& str, char separator) {
    int startIndex = 0, endIndex = 0;
    std::vector<std::string> strings;

    for (size_t i = 0; i <= str.size(); i++) {
      // If we reached the end of the word or the end of the input.
      if (str[i] == separator || i == str.size()) {
        endIndex = i;
        std::string temp;
        temp.append(str, startIndex, endIndex - startIndex);
        strings.push_back(temp);
        startIndex = endIndex + 1;
      }
    }
    return strings;
  }

  std::optional<std::filesystem::path> find_file_in_dir(const std::filesystem::path& base_path,
                                                        std::string_view filename) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(base_path)) {
      if (entry.is_directory()) continue;
      if (entry.path().filename().string() == filename) return entry.path();
      if (entry.path().stem().string() == filename) return entry.path();
    }
    return std::nullopt;
  };

  std::filesystem::path get_file_path(std::string_view filename) {
    auto filepath = find_file_in_dir(get_base_path() + "/data", filename);
    if (!filepath.has_value()) {
      std::string msg = "File not found: " + std::string(filename);
      throw std::runtime_error(msg);
    }
    return filepath.value();
  }

  TLE TLE::FromLines(const std::string& line1, const std::string& line2, const std::string& line3) {
    TLE tle;
    tle.name = trim_right_copy(line1);
    if (line1.substr(0, 3) == "GPS") {
      tle.prn = stod(split_string(line1, '(')[1].substr(4, 2));
    } else if (line1.substr(0, 3) == "BEI") {
      std::vector<std::string> split = split_string(line1, '(');
      if (split[1].substr(0, 1) == "C") {
        tle.prn = stod(split[1].substr(1, 2));
      } else {
        tle.prn = 0;
      }
    } else if ((line1.substr(0, 3) == "GSA") && (line1.substr(4, 2) != "01")) {
      std::string tmp = split_string(line1, '(')[1];
      std::string tmp2 = split_string(tmp, ' ')[1];
      std::string tmp3 = split_string(tmp2, ')')[0];
      if (tmp3.size() == 1) {  // ( ex. 8) )
        tle.prn = stod(tmp3.substr(0, 1));
      } else if (tmp3.size() == 2) {
        tle.prn = stod(tmp3.substr(0, 2));
      }
    } else if (line1.substr(0, 3) == "COS") {
      tle.prn = stod(split_string(line1, '(')[1].substr(0, 3));
    } else if (line1.substr(0, 3) == "QZS") {
      tle.prn = stod(line1.substr(4, 1));  // 1-4

    } else {
      tle.prn = -1;
    }
    tle.epoch_year = stod(line2.substr(18, 2));
    tle.epoch_day = stod(line2.substr(20, 12));
    tle.bstar = stod(line2.substr(53, 8));
    tle.inclination = stod(line3.substr(8, 8));
    tle.raan = stod(line3.substr(17, 8));
    tle.eccentricity = stod("0." + line3.substr(26, 7));
    tle.arg_perigee = stod(line3.substr(34, 8));
    tle.mean_anomaly = stod(line3.substr(43, 8));
    tle.mean_motion = stod(line3.substr(52, 11));

    // compute TAI from epoch
    // std::string fullyear_string = "20" + line2.substr(18, 2);
    // Real epoch_year_start_tai = spice::String2TAI(fullyear_string + "/01/01
    // 00:00:00 UTC");
    double mjd = gregorian_to_mjd(2000 + tle.epoch_year, 1, 1, 0, 0, 0);
    double epoch_year_start_utc = mjd_to_tj2000(mjd);
    double epoch_utc = epoch_year_start_utc + tle.epoch_day * SECS_DAY;
    tle.epoch_utc = epoch_utc;
    return tle;
  };

  std::vector<TLE> TLE::FromFile(const std::string_view filename) {
    std::filesystem::path path = get_file_path(filename);
    std::ifstream input_file(path);
    if (!input_file.is_open()) {
      throw std::runtime_error("Could not open file " + std::string(filename));
    }
    std::vector<TLE> tles;
    std::string line1, line2, line3;
    while (getline(input_file, line1) && getline(input_file, line2) && getline(input_file, line3)) {
      TLE tle = TLE::FromLines(line1, line2, line3);
      if (tle.prn != -1) {
        tles.push_back(tle);
      }
    };
    return tles;
  };

}  // namespace pecsim
