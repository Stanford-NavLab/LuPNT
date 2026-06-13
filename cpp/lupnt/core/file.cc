#include "lupnt/core/file.h"

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"

namespace lupnt {
  std::filesystem::path GetBaseDir() {
    const char* base_path_env = std::getenv("LUPNT_PATH");
    LUPNT_CHECK(base_path_env, "Environment variable LUPNT_PATH is not set.", "File");
    return std::filesystem::path(base_path_env);
  }

  std::filesystem::path GetDataPath() {
    const char* data_path_env = std::getenv("LUPNT_DATA_PATH");
    if (data_path_env == nullptr) {
      std::string msg = "Environment variable LUPNT_DATA_PATH is not set.";
      LUPNT_CHECK(false, msg, "File");
    }
    return std::filesystem::path(data_path_env);
  }

  std::filesystem::path GetOutputDir(const std::string& dirname) {
    const char* output_path_env = std::getenv("LUPNT_OUTPUT_PATH");
    std::filesystem::path output_path;
    if (output_path_env == nullptr) {
      output_path = GetDataPath() / "output";
    } else {
      output_path = std::filesystem::path(output_path_env);
    }
    if (!dirname.empty()) output_path /= dirname;
    if (!std::filesystem::exists(output_path)) std::filesystem::create_directories(output_path);
    return output_path;
  }

  std::filesystem::path SetOutputDir(const std::string& dirname) {
    auto output_dir = GetOutputDir(dirname);
    Logger::Info(fmt::format("Output directory: {}", output_dir.string()), "File", 0.0);
    return output_dir;
  }

  std::filesystem::path GetInputPath(const std::string& input_file) {
    auto input_path_env = GetDataPath() / "..";
    return FindFileInDir(input_path_env, input_file).value();
  }

  YAML::Node GetYamlConfig(const std::filesystem::path& filepath) {
    return YAML::LoadFile(filepath);
  }

  H5Easy::File GetH5File(std::filesystem::path cache_path, bool recompute) {
    if (recompute && std::filesystem::exists(cache_path)) std::filesystem::remove(cache_path);
    auto open_mode = (recompute || !std::filesystem::exists(cache_path))
                         ? H5Easy::File::OpenOrCreate
                         : H5Easy::File::ReadWrite;
    H5Easy::File cache_file(cache_path, open_mode);
    return cache_file;
  }

  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                     std::string_view filename) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(base_path)) {
      if (entry.is_directory()) continue;
      if (entry.path().filename().string() == filename) return entry.path();
      if (entry.path().stem().string() == filename) return entry.path();
    }
    return std::nullopt;
  };

  std::filesystem::path GetFilePath(std::string_view filename) {
    auto filepath = FindFileInDir(GetDataPath(), filename);
    if (!filepath.has_value()) {
      std::string msg = "File not found: " + std::string(filename);
      throw std::runtime_error(msg);
    }
    return filepath.value();
  }

  std::filesystem::path GetCspiceKernelDir() { return GetDataPath() / "ephemeris"; }
  std::filesystem::path GetAsciiKernelDir() { return GetDataPath() / "ephemeris" / "ascii"; }
  std::filesystem::path GetConfigFileDir() { return GetDataPath() / ".." / "config"; }

  std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime() {
    return std::chrono::high_resolution_clock::now();
  }

  std::string PrintDuration(const std::chrono::duration<double>& duration, int precision) {
    auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    auto hours = total_seconds / 3600;
    auto minutes = (total_seconds % 3600) / 60;
    auto seconds = total_seconds % 60;
    double fractional_seconds = duration.count() - total_seconds;  // Get decimal part

    std::ostringstream oss;
    if (hours > 0) oss << hours << "h ";
    if (minutes > 0) oss << minutes << "m ";

    // Format seconds with precision
    oss << std::fixed << std::setprecision(precision) << (seconds + fractional_seconds) << "s";

    return oss.str();
  }

  template <typename T> T OpenFile(const std::filesystem::path& filepath) {
    T file(filepath);
    LUPNT_CHECK(file.is_open(), "Unable to open file: " + filepath.string(), "File");
    return file;
  }
  template std::ifstream OpenFile<std::ifstream>(const std::filesystem::path& filepath);
  template std::ofstream OpenFile<std::ofstream>(const std::filesystem::path& filepath);

  size_t CountLines(const std::filesystem::path& filepath) {
    std::ifstream file = OpenFile<std::ifstream>(filepath);
    size_t line_count = 0;
    std::string line;
    while (std::getline(file, line)) ++line_count;
    file.close();
    return line_count;
  }

  void Dump(H5Easy::File& file, const std::string& name, Real value, H5Easy::DumpMode mode) {
    double tmp = value.val();
    H5Easy::dump(file, name, tmp, mode);
  }

}  // namespace lupnt
