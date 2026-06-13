/**
 * @file tec/klobuchar.h
 * @author Keidai Iiyama
 * @brief This file contains the interface for the Klobuchar model.
 * @version 0.1
 * @date 2025-02-17
 */

#pragma once

#include <cmath>

#include "lupnt/environment/plasma/core/definitions.h"

namespace pecsim {

  double Klobucher(double t_gps, double elevation, double azimuth, double latitude_u,
                   double longitude_u, double freq_Hz, const Vec4d& alpha, const Vec4d& beta);

}  // namespace pecsim
