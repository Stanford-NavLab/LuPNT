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
  Real x_pole;
  Real y_pole;
  Real UT1_UTC;
  Real LOD;
  Real dPsi;
  Real dEps;
  Real dx_pole;
  Real dy_pole;
  Real TAI_UTC;
};

// Function to count lines in the file
size_t CountLines(const std::string& filename);

// Function to load EOP data from the file
std::shared_ptr<EOPData> LoadEOPData(const std::string& filename);

// Function to manage IERS time and polar motion data
EOPResult InterpolateEOPData(const std::shared_ptr<EOPData>& eop_data,
                             Real mjdUTC, bool interpolate);

}  // namespace lupnt