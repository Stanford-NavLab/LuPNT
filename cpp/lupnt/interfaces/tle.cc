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

#include "lupnt/interfaces/tle.h"

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <memory>
#include <string>

#include "lupnt/agents/agent.h"
#include "lupnt/agents/satellite.h"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"
#include "lupnt/core/string_utils.h"
#include "lupnt/states/state.h"

namespace lupnt {

  TLE TLE::FromLines(const std::string& line1, const std::string& line2, const std::string& line3) {
    auto parse_bstar = [](const std::string& s) {
      std::string mant = s.substr(0, 5);
      std::string exp = s.substr(5, 3);
      return std::stod(mant + "E" + exp);
    };
    auto parse_year = [](double yy) { return (yy < 57 ? 2000 : 1900) + static_cast<int>(yy); };

    TLE tle;
    tle.name = boost::trim_right_copy(line1);
    tle.epoch_year = parse_year(stod(line2.substr(18, 2)));
    tle.epoch_day = stod(line2.substr(20, 12));
    tle.bstar = parse_bstar(line2.substr(53, 8));
    tle.inclination = stod(line3.substr(8, 8));
    tle.raan = stod(line3.substr(17, 8));
    tle.eccentricity = stod("0." + line3.substr(26, 7));
    tle.arg_perigee = stod(line3.substr(34, 8));
    tle.mean_anomaly = stod(line3.substr(43, 8));
    tle.mean_motion = stod(line3.substr(52, 11));

    Real epoch_year_start_utc = GregorianToTime(tle.epoch_year, 1, 1, 0, 0, 0);
    Real epoch_year_start_tai = UtcToTai(epoch_year_start_utc);

    double epoch_tai = epoch_year_start_tai.val() + tle.epoch_day * SECS_DAY;
    tle.epoch_tai = epoch_tai;
    tle.prn = GetPrn(tle.name);
    return tle;
  };

  std::vector<TLE> TLE::FromFile(const std::string_view filename) {
    std::filesystem::path path = GetFilePath(filename);
    std::ifstream input_file(path);
    LUPNT_CHECK(input_file.is_open(), "Could not open file " + std::string(filename), "TLE");
    std::vector<TLE> tles;
    std::string line1, line2, line3;
    while (getline(input_file, line1) && getline(input_file, line2) && getline(input_file, line3)) {
      TLE tle = TLE::FromLines(line1, line2, line3);
      if (tle.prn != -1) tles.push_back(tle);
    };
    return tles;
  };

  std::vector<Ptr<Agent>> LoadTleFile(const std::string& filename) {
    std::vector<Ptr<Agent>> satellites;

    for (auto tle : TLE::FromFile(filename)) {
      // Classical orbital elements
      Real T = SECS_DAY / tle.mean_motion;
      Real a = pow((T * T * GM_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
      Real e = tle.eccentricity;
      Real i = tle.inclination * RAD;
      Real Omega = tle.raan * RAD;
      Real w = tle.arg_perigee * RAD;
      Real rad_per_sec = tle.mean_motion * 2.0 * PI / SECS_DAY;
      Real dt = GetLupntEpoch() - tle.epoch_tai;
      Real M = WrapToPi(tle.mean_anomaly * RAD + rad_per_sec * dt);

      // Convert to Cartesian: a Satellite carries a Cartesian state (its
      // dynamics integrate rv), so the classical elements must be converted
      // before SetState -- passing the ClassicalOE directly leaves the agent
      // with a non-Cartesian state that breaks downstream propagation.
      State coe = ClassicalOE({a, e, i, Omega, w, M}, Frame::GCRF);
      State rv = ClassicalToCart(coe, GM_EARTH);

      // Create the spacecraft
      auto sat = MakePtr<Satellite>();
      sat->SetTime(0.0);
      sat->SetState(rv);
      sat->SetName(tle.name);
      satellites.push_back(std::move(sat));
    }

    Logger::Debug(
        fmt::format("Loaded {} satellites from TLE file: {}", satellites.size(), filename), "TLE",
        0.0);
    return satellites;
  }

  int GetPrn(const std::string& name) {
    // PRN
    if (name.size() < 3) return -1;

    std::string prefix = name.substr(0, 3);
    std::vector<std::string> split;

    if (prefix == "GPS") {
      // GPS BIIR-2  (PRN 13)  -> 13
      return std::stoi(SplitString(name, '(')[1].substr(4, 2));
    }
    if (prefix == "BEI") {
      // BEIDOU-2 M3 (C11)     -> 11
      split = SplitString(name, '(');
      return (split[1][0] == 'C') ? std::stoi(split[1].substr(1, 2)) : 0;
    }
    if (prefix == "GSA") {
      // GSAT0202 (GALILEO 6)  -> 6
      if (name.substr(4, 2) == "01") return -1;
      split = SplitString(SplitString(SplitString(name, '(')[1], ' ')[1], ')');
      return split[0].size() <= 2 ? std::stoi(split[0]) : -1;
    }
    if (prefix == "COS") {
      // COSMOS 2501 (702K)    -> 702
      return std::stoi(SplitString(name, '(')[1].substr(0, 3));
    }
    if (prefix == "QZS") {
      // QZS-1R (QZSS/PRN 196) -> 196
      return std::stoi(name.substr(4, 1));
    }

    std::cout << prefix << std::endl;
    std::cout << name.substr(4, 2) << std::endl;
    std::cout << (prefix == "GSA") << std::endl;
    std::cout << (name.substr(4, 2) == "01") << std::endl;
    LUPNT_CHECK(false, "Unknown GNSS PRN format for name: " + name, "TLE");
    return -1;
  }

}  // namespace lupnt
