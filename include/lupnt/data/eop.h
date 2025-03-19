#pragma once

#include <lupnt/core/constants.h>

#include <memory>
#include <mutex>
#include <string>

namespace lupnt {

  // Struct to hold EOP data
  struct EopFileData {
    VecXi years;
    VecXi months;
    VecXi days;
    VecXd mjds_utc;
    VecXd x;
    VecXd y;
    VecXd ut1_utc;
    VecXd lod;
    VecXd dpsi;
    VecXd deps;
    VecXd xErr;
    VecXd yErr;
    VecXd ut1_utc_err;
    VecXd lod_err;
    VecXd dpsi_err;
    VecXd deps_err;
  };

  struct EopData {
    Real x_pole;
    Real y_pole;
    Real ut1_utc;
    Real lod;
    Real dpsi;
    Real deps;
    Real dx_pole;
    Real dy_pole;
    Real tai_utc;
  };

  void LoadEopFileData(const std::filesystem::path& filepath);

  // Function to manage IERS time and polar motion data
  EopData GetEopData(Real mjd_utc);

  Real GetUt1UtcDifference(Real mjd_utc);

}  // namespace lupnt
