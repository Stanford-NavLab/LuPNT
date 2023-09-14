

#pragma once

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>

namespace LPT {
// set here the absolute path to your lunanet simulator
static const std::filesystem::path BASEPATH("/path/to/this/LuPNT");

static std::optional<std::filesystem::path> FindFileInDir(
    const std::filesystem::path& basePath, const std::string& filename) {
  for (const auto& entry :
       std::filesystem::recursive_directory_iterator(basePath)) {
    if (entry.path().stem().string() == filename) {
      return entry.path();
    }
  }
  return std::nullopt;
};

}  // namespace LPT