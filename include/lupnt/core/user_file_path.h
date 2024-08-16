#pragma once

#include <assert.h>

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>

namespace lupnt {

  std::filesystem::path GetDataPath();
  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& basePath,
                                                     const std::string& filename);
  std::filesystem::path GetFilePath(const std::string& filename);
  std::filesystem::path GetCspiceKernelDir();
  std::filesystem::path GetAsciiKernelDir();

}  // namespace lupnt
