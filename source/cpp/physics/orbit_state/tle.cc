/**
 * @file TLE.cpp
 * @author Stanford NAV LAB
 * @brief  Routines to process TLE files
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/physics/orbit_state/tle.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "lupnt/numerics/string_utils.h"
#include "lupnt/physics/orbit_state/orbit_states.h"
#include "lupnt/physics/spice_interface.h"

namespace lupnt {
  TLE TLE::FromLines(const std::string& line1, const std::string& line2, const std::string& line3) {
    TLE tle;
    if (line1.substr(0, 3) == "GPS") {
      tle.name = "GPS";
      tle.prn = stod(SplitString(line1, '(')[1].substr(4, 2));
    } else if (line1.substr(0, 3) == "BEI") {
      tle.name = "BEIDOU";
      std::vector<std::string> split = SplitString(line1, '(');
      if (split[1].substr(0, 1) == "C") {
        tle.prn = stod(split[1].substr(1, 2));
      } else {
        tle.prn = 0;
      }
    } else if (line1.substr(0, 3) == "GSA") {
      tle.name = "GALILEO";
      tle.prn = stod(SplitString(line1, '(')[1].substr(5, 2));
    } else if (line1.substr(0, 3) == "COS") {
      tle.name = "GLONASS";
      tle.prn = stod(SplitString(line1, '(')[1].substr(0, 3));
    } else {
      tle.name = line1;
      tle.prn = 0;
    }
    tle.epochYear = stod(line2.substr(18, 2));
    tle.epochDay = stod(line2.substr(20, 12));
    tle.bstar = stod(line2.substr(53, 8));
    tle.inclination = stod(line3.substr(8, 8));
    tle.raan = stod(line3.substr(17, 8));
    tle.eccentricity = stod("0." + line3.substr(26, 7));
    tle.argPerigee = stod(line3.substr(34, 8));
    tle.meanAnomaly = stod(line3.substr(43, 8));
    tle.meanMotion = stod(line3.substr(52, 11));

    // compute TAI from epoch
    std::string fullyear_string = "20" + line2.substr(18, 2);
    Real epochYearStartTAI = spice::String2TAI(fullyear_string + "/01/01 00:00:00 UTC");
    double epochTAI = epochYearStartTAI.val() + tle.epochDay * SECS_DAY;
    tle.epochTAI = epochTAI;

    return tle;
  };

  std::vector<TLE> TLE::FromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
      throw std::runtime_error("Could not open file " + filename);
    }
    std::vector<TLE> tles;
    std::string line1, line2, line3;
    while (getline(inputFile, line1) && getline(inputFile, line2) && getline(inputFile, line3)) {
      TLE tle = TLE::FromLines(line1, line2, line3);
      tles.push_back(tle);
    };
    return tles;
  };
}  // namespace lupnt
