#pragma once
#include <memory>
#include <string>
#include <vector>

namespace lupnt {
  class TLE {
  public:
    std::string name;
    double epoch_year;    // [yr] Last two digits of year
    double epoch_day;     // [day] Day of the year and fractional portion of the day
    double epoch_tai;     // [s TAI]
    double bstar;         // [1/R_EARTH] B*, the drag term, or radiation pressure coefficient
    double inclination;   // [deg] Inclination
    double raan;          // [deg] Right ascension of the ascending node
    double eccentricity;  // [-] Eccentricity
    double arg_perigee;   // [deg] Argument of perigee
    double mean_anomaly;  // [deg] Mean anomaly
    double mean_motion;   // [revs/day] Mean motion
    int prn;              // PRN number
    static TLE FromLines(const std::string& line1, const std::string& line2,
                         const std::string& line3);
    static std::vector<TLE> FromFile(const std::string& filename);
  };

}  // namespace lupnt
