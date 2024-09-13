#pragma once
#include <memory>
#include <string>
#include <vector>

namespace lupnt {
  class TLE {
  public:
    std::string name;
    double epoch_year;
    double epoch_day;
    double epoch_tai;
    double bstar;
    double inclination;
    double raan;
    double eccentricity;
    double arg_perigee;
    double mean_anomaly;
    double mean_motion;
    int prn;
    static TLE FromLines(const std::string& line1, const std::string& line2,
                         const std::string& line3);
    static std::vector<TLE> FromFile(const std::string& filename);
  };

}  // namespace lupnt
