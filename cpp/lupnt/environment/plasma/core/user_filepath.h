#pragma once

#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#if defined(__APPLE__)
#  include <mach-o/dyld.h>
#elif defined(__linux__)
#  include <unistd.h>
#endif

namespace pecsim {

  inline bool has_runtime_data(const std::filesystem::path& base_path) {
    std::error_code ec;
    return std::filesystem::exists(base_path / "data" / "iri", ec);
  }

  inline std::string executable_path() {
#if defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::vector<char> buffer(size);
    if (_NSGetExecutablePath(buffer.data(), &size) == 0) {
      return std::filesystem::weakly_canonical(buffer.data()).string();
    }
#elif defined(__linux__)
    std::vector<char> buffer(4096, '\0');
    ssize_t size = readlink("/proc/self/exe", buffer.data(), buffer.size() - 1);
    if (size > 0) {
      buffer[static_cast<size_t>(size)] = '\0';
      return std::string(buffer.data());
    }
#endif
    return std::string();
  }

  inline std::string infer_base_path_from_executable() {
    std::filesystem::path current = executable_path();
    if (current.empty()) {
      return std::string();
    }

    current = current.parent_path();
    for (int i = 0; i < 8 && !current.empty(); ++i) {
      if (has_runtime_data(current)) {
        return current.string();
      }
      current = current.parent_path();
    }

    return std::string();
  }

  inline std::string& mutable_base_path() {
    static std::string base_path = []() {
      const char* env = std::getenv("PECSIMPY_BASE_PATH");
      if (env != nullptr) {
        return std::string(env);
      }
      return infer_base_path_from_executable();
    }();
    return base_path;
  }

  inline void set_base_path(const std::string& path) { mutable_base_path() = path; }

  inline const std::string& get_base_path() {
    const char* env = std::getenv("PECSIMPY_BASE_PATH");
    if (env != nullptr && mutable_base_path() != env) {
      mutable_base_path() = env;
    }
    return mutable_base_path();
  }

}  // namespace pecsim
