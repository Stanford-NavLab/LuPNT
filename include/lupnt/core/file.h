/**
 * @file file.h
 * @author Stanford NAV LAB
 * @brief File access utils
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <filesystem>
#include <fstream>
#include <highfive/H5Easy.hpp>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {
  using H5Easy::DataSet;
  using H5Easy::dump;
  using H5Easy::DumpMode;
  using H5Easy::DumpOptions;
  using H5Easy::File;
  using H5Easy::load;

  size_t CountLines(const std::filesystem::path& filepath);

}  // namespace lupnt
