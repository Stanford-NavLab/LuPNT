#pragma once
#include <memory>
#include <string>
#include <vector>
namespace lupnt {
  class TLE {
  public:
    std::string name;
    double epochYear;
    double epochDay;
    double epochTAI;
    double bstar;
    double inclination;
    double raan;
    double eccentricity;
    double argPerigee;
    double meanAnomaly;
    double meanMotion;
    int prn;

    static TLE FromLines(const std::string &line1, const std::string &line2,
                         const std::string &line3);
    static std::vector<TLE> FromFile(const std::string &filename);
  };

}  // namespace lupnt
