#pragma once

#include <assert.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>

namespace lupnt {

  std::filesystem::path GetDataPath();
  std::filesystem::path GetOutputPath(std::string output_dir);
  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                     std::string_view filename);
  std::filesystem::path GetFilePath(std::string_view filename);
  std::filesystem::path GetCspiceKernelDir();
  std::filesystem::path GetAsciiKernelDir();

  std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime();
  std::string PrintDuration(const std::chrono::duration<double>& duration);

}  // namespace lupnt
