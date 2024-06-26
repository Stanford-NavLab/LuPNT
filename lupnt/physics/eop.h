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
  VectorXi years;
  VectorXi months;
  VectorXi days;
  VectorXi mjds;
  VectorXd x;
  VectorXd y;
  VectorXd ut1_utc;
  VectorXd lod;
  VectorXd dPsi;
  VectorXd dEps;
  VectorXd xErr;
  VectorXd yErr;
  VectorXd ut1_utcErr;
  VectorXd lodErr;
  VectorXd dPsiErr;
  VectorXd dEpsErr;
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