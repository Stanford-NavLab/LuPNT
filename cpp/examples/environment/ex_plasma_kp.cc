/**
 * @file test_kp.cpp
 * @author Keidai Iiyama
 * @brief Test Kp index retrieval
 * @version 0.1
 * @date 2025-02-17
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <iostream>

#include "lupnt/environment/plasma/plasma.h"

using namespace std;
using namespace pecsim;

int main() {
  std::vector<int> years = {2000, 2005, 2010, 2015, 2020, 2024};
  std::vector<int> doys = {1, 100, 200, 250, 300, 365};
  std::vector<int> hours = {0, 2, 8, 11, 18, 23};
  int n = years.size();

  DateTime datetime;
  datetime.min = 25;
  datetime.sec = 39.1;

  for (int i = 0; i < n; i++) {
    datetime.year = years[i];
    datetime.doy = doys[i];
    datetime.hour = hours[i];
    double kp = get_kp_index(datetime);
    cout << "Year: " << datetime.year << " Day of Year: " << datetime.doy
         << " Hour: " << datetime.hour << " Kp: " << kp << endl;
  }

  return 0;
}
