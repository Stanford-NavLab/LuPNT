#pragma once

#include <lupnt/core/constants.h>

#include <memory>
#include <mutex>
#include <string>

namespace lupnt {

  // Struct to hold EOP data
  struct IauSofaFileData {
    VecXd jd_tt;
    VecXd X;
    VecXd Y;
    VecXd s;
  };

  struct IauSofaData {
    Real X;
    Real Y;
    Real s;
  };

  void LoadIauSofaFileData(const std::filesystem::path& filepath);

  IauSofaData GetIauSofaData(Real jd_tt);

}  // namespace lupnt
