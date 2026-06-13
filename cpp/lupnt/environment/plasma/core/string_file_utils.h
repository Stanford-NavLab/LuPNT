/**
 * @file string_utils.h
 * @author Keidai Iiyama
 * @brief This file contains utility functions for string and file manipulation.
 * @version 0.1
 * @date 2025-02-17
 */

#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace pecsim {
  std::string trim_right_copy(const std::string& input);
  std::vector<std::string> split_string(const std::string& str, char separator);
  std::filesystem::path get_file_path(const std::string& filename);
  std::optional<std::filesystem::path> find_file_in_dir(const std::filesystem::path& base_path,
                                                        std::string_view filename);

  /**
   * @brief Class representing a Two-Line Element (TLE) set.
   *
   * This class provides methods to parse TLE data from strings or files.
   */
  class TLE {
  public:
    std::string name;
    double epoch_year;    // [yr] Last two digits of year
    double epoch_day;     // [day] Day of the year and fractional portion of the day
    double epoch_utc;     // [s TAI]
    double bstar;         // [1/R_EARTH] B*, the drag term, or radiation pressure
                          // coefficient
    double inclination;   // [deg] Inclination
    double raan;          // [deg] Right ascension of the ascending node
    double eccentricity;  // [-] Eccentricity
    double arg_perigee;   // [deg] Argument of perigee
    double mean_anomaly;  // [deg] Mean anomaly
    double mean_motion;   // [revs/day] Mean motion
    int prn;              // PRN number
    static TLE FromLines(const std::string& line1, const std::string& line2,
                         const std::string& line3);
    static std::vector<TLE> FromFile(const std::string_view filename);
  };

}  // namespace pecsim
