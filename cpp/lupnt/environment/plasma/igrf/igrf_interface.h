/**
 * @file igrf_interface.h
 * @author Keidai Iiyama
 * @brief This file contains the interface for the IGRF model.
 * @version 0.1
 * @date 2025-02-06
 */

#pragma once

#include <vector>

namespace pecsim {

  std::vector<double> igrf14(double lat_deg, double lon_deg, double r_km, double decimal_year);

}
