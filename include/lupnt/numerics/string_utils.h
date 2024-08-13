/**
 * @file string_utils.h
 * @author Stanford NAV LAB
 * @brief String util functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace lupnt {

  std::vector<std::string> SplitString(const std::string& str, char separator);
  std::vector<std::vector<std::string>> ReadCSV(std::filesystem::path fname);

}  // namespace lupnt
