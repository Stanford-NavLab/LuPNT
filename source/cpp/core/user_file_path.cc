#include "lupnt/core/user_file_path.h"

namespace lupnt {
  std::filesystem::path GetDataPath() {
    const char* dataPathEnv = std::getenv("LUPNT_DATA_PATH");
    assert(dataPathEnv != nullptr && "Environment variable LUPNT_DATA_PATH is not set.");
    return std::filesystem::path(dataPathEnv);
  }

  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& basePath,
                                                     const std::string& filename) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(basePath)) {
      if (entry.path().filename().string() == filename) {
        return entry.path();
      }
    }
    return std::nullopt;
  };

  std::filesystem::path GetFilePath(const std::string& filename) {
    auto filepath = FindFileInDir(GetDataPath(), filename);
    assert(filepath.has_value() && "File not found.");
    return filepath.value();
  }

  std::filesystem::path GetCspiceKernelDir() { return GetDataPath() / "ephemeris"; }
  std::filesystem::path GetAsciiKernelDir() { return GetDataPath() / "ephemeris" / "ascii"; }

}  // namespace lupnt
