#pragma once

#include <lupnt/core/constants.h>

#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <memory>
#include <mutex>
#include <string>

namespace lupnt {

// Struct to hold EOP data
struct EOPData {
  VecXi years;
  VecXi months;
  VecXi days;
  VecXi mjds;
  VecXd x;
  VecXd y;
  VecXd ut1_utc;
  VecXd lod;
  VecXd dPsi;
  VecXd dEps;
  VecXd xErr;
  VecXd yErr;
  VecXd ut1_utcErr;
  VecXd lodErr;
  VecXd dPsiErr;
  VecXd dEpsErr;
};

struct EOPResult {
  real x_pole;
  real y_pole;
  real UT1_UTC;
  real LOD;
  real dPsi;
  real dEps;
  real dx_pole;
  real dy_pole;
  real TAI_UTC;
};

// Function to count lines in the file
size_t CountLines(const std::string& filename);

// Function to load EOP data from the file
std::shared_ptr<EOPData> LoadEOPData(const std::string& filename);

// Function to manage IERS time and polar motion data
EOPResult InterpolateEOPData(const std::shared_ptr<EOPData>& eop_data,
                             real mjdUTC, bool interpolate);

}  // namespace lupnt