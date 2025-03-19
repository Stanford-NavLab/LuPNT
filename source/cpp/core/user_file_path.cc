#include "lupnt/core/user_file_path.h"

namespace lupnt {
  std::filesystem::path GetDataPath() {
    const char* data_path_env = std::getenv("LUPNT_DATA_PATH");
    assert(data_path_env != nullptr && "Environment variable LUPNT_DATA_PATH is not set.");
    return std::filesystem::path(data_path_env);
  }

  std::filesystem::path GetOutputPath(std::string output_dir) {
    const char* output_path_env = std::getenv("LUPNT_OUTPUT_PATH");
    std::filesystem::path output_path;
    if (output_path_env == nullptr) {
      output_path = GetDataPath() / "output";
    } else {
      output_path = std::filesystem::path(output_path_env);
    }
    if (!output_dir.empty()) {
      output_path /= output_dir;
    }

    if (!std::filesystem::exists(output_path)) {
      std::filesystem::create_directories(output_path);
    }
    return output_path;
  }

  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                     std::string_view filename) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(base_path)) {
      if (entry.path().filename().string() == filename) {
        return entry.path();
      }
    }
    return std::nullopt;
  };

  std::filesystem::path GetFilePath(std::string_view filename) {
    auto filepath = FindFileInDir(GetDataPath(), filename);
    assert(filepath.has_value() && "File not found.");
    return filepath.value();
  }

  std::filesystem::path GetCspiceKernelDir() { return GetDataPath() / "ephemeris"; }
  std::filesystem::path GetAsciiKernelDir() { return GetDataPath() / "ephemeris" / "ascii"; }

  std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime() {
    return std::chrono::high_resolution_clock::now();
  }
  // print duration between start_time and end_time
  std::string PrintDuration(const std::chrono::duration<double>& duration) {
    auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    auto hours = total_seconds / 3600;
    auto minutes = (total_seconds % 3600) / 60;
    auto seconds = total_seconds % 60;

    std::string str;
    if (hours > 0) str += std::to_string(hours) + "h ";
    if (minutes > 0) str += std::to_string(minutes) + "m ";
    str += std::to_string(seconds) + "s ";
    return str;
  }
}  // namespace lupnt
