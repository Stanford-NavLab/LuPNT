/**
 * @file igrf_interface.cpp
 * @author Keidai Iiyama
 * @brief This file contains the implementation of the IGRF model interface.
 * @version 0.1
 * @date 2025-02-06
 */

#include "lupnt/environment/plasma/igrf/igrf_interface.h"

namespace pecsim {
  extern "C" {
  void igrf14syn_(int* isv, double* date, int* itype, double* alt, double* colat, double* elong,
                  double* x, double* y, double* z, double* f);
  }

  std::vector<double> igrf14(double lat_rad, double lon_rad, double r_km, double decimal_year) {
    int isv = 0;    // 0 = magnetic field, not secular variation
    int itype = 2;  // 1 = geodetic, 2 = geocentric
    double x, y, z, f;

    double lat_deg = lat_rad * 180.0 / M_PI;  // Convert latitude from radians to degrees
    double lon_deg = lon_rad * 180.0 / M_PI;  // Convert longitude from radians to degrees

    double colat_deg = 90.0 - lat_deg;  // Convert latitude to colatitude
    double elong_deg = lon_deg;         // Longitude remains the same
    if (elong_deg < 0) {
      elong_deg += 360.0;  // Ensure longitude is in the range [0, 360)
    }

    // Note: Fortran passes everything by reference
    igrf14syn_(&isv, &decimal_year, &itype, &r_km, &colat_deg, &elong_deg, &x, &y, &z, &f);

    std::vector<double> result(3);
    result[0] = x;  // X component of the magnetic field in nT (north)
    result[1] = y;  // Y component of the magnetic field in nT (east)
    result[2] = z;  // Z component of the magnetic field in nT (vertical)

    return result;
  }

}  // namespace pecsim
