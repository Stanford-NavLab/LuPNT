/**
 *  @file kp.cpp
 *  @author Keidai Iiyama
 *  @brief This file contains the implementation of the Kp index retrieval.
 *  @version 0.1
 *  @date 2025-02-17
 */

#include "lupnt/environment/plasma/env/kp.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "lupnt/environment/plasma/core/user_filepath.h"

namespace pecsim {

  double get_kp_index(DateTime datetime) {
    // get the Kp index from the data
    std::string csv_dir = get_base_path() + "/data/kp/csv";
    double kp = 0.0;

    int year = datetime.year;
    int doy = datetime.doy;
    int hour = datetime.hour;
    int min = datetime.min;
    double sec = datetime.sec;

    double hourd = hour + min / 60.0 + sec / 3600.0;

    int mm, dd;
    doy_to_mmdd(year, doy, mm, dd);

    // closest hour from [0, 3, 6, 9, 12, 15, 18, 21]
    int closest_hour = 3 * std::round(hourd / 3.0);

    // mod24
    closest_hour = closest_hour % 24;

    // Read the Kp index from the CSV file
    // The CSV file should have the following columns:
    // year, mm, dd, hh, kp
    // The CSV file should be named as kp_YYYY.csv
    // where YYYY is the year

    std::string filename = csv_dir + "/kp_" + std::to_string(year) + ".csv";

    // Read the CSV file
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::runtime_error("Could not open the Kp index file.");
    }

    std::string line;
    std::vector<std::string> row;

    // get line
    while (std::getline(file, line)) {
      row.clear();
      std::string cell;
      std::stringstream lineStream(line);

      while (std::getline(lineStream, cell, ',')) {
        row.push_back(cell);
      }

      // ignore line 1
      if (row[0] == "yyyy") {
        continue;
      }
      if (std::stoi(row[0]) == year && std::stoi(row[1]) == mm && std::stoi(row[2]) == dd
          && std::stoi(row[3]) == closest_hour) {
        kp = std::stod(row[4]);
        break;
      }
    }

    return kp;
  }

}  // namespace pecsim
