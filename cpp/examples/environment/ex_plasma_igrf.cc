/**
 * @file ex_plasma_igrf.cc
 * @author Keidai Iiyama
 * @brief Evaluate the IGRF14 magnetic field model at one geodetic point.
 * @version 0.1
 * @date 2025-02-15
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;
using namespace std;

int main() {
  // IGRF expects geocentric radius in km and latitude/longitude in radians.
  double lat_deg = 45.0;         // Latitude in degrees
  double lon_deg = -80.0;        // Longitude in degrees
  double r_km = 2 * 6371.0;      // Radius in kilometers (Earth's radius)
  double decimal_year = 2024.5;  // Decimal year for the IGRF model

  double lat_rad = lat_deg * M_PI / 180.0;  // Convert latitude to radians
  double lon_rad = lon_deg * M_PI / 180.0;  // Convert longitude to radians

  std::vector<double> Bxyz = igrf14(lat_rad, lon_rad, r_km, decimal_year);

  cout << "IGRF Magnetic Field Components at (" << lat_deg << "°, " << lon_deg << "°, " << r_km
       << " km):" << endl;
  cout << "  X: " << Bxyz[0] << " nT (North)" << endl;
  cout << "  Y: " << Bxyz[1] << " nT (East)" << endl;
  cout << "  Z: " << Bxyz[2] << " nT (Vertical)" << endl;
}
