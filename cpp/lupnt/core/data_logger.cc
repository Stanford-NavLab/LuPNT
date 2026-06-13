#include "lupnt/core/data_logger.h"

#include <fmt/format.h>

#include <highfive/H5Easy.hpp>

#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"

namespace lupnt {

  // Use function-local static for thread-safe initialization
  DataLogger::DoubleMap& GetDoubleLogs() {
    static DataLogger::DoubleMap double_logs;
    return double_logs;
  }

  DataLogger::StringMap& GetStringLogs() {
    static DataLogger::StringMap string_logs;
    return string_logs;
  }

  DataLogger::MatrixMap& GetMatrixLogs() {
    static DataLogger::MatrixMap matrix_logs;
    return matrix_logs;
  }

  std::mutex& GetMutex() {
    static std::mutex mutex;
    return mutex;
  }

  std::unique_ptr<H5Easy::File>& GetFile() {
    static std::unique_ptr<H5Easy::File> file = nullptr;
    return file;
  }

  // Set output file for logging
  void DataLogger::SetOutputFile(const std::filesystem::path& filename) {
    Logger::Info("Output file: {}", filename.string());
    std::lock_guard<std::mutex> lock(GetMutex());

    if (std::filesystem::exists(filename)) {
      Logger::Debug(fmt::format("Removing existing output file: {}", filename.string()),
                    "DataLogger");
      std::filesystem::remove(filename);
    }

    GetFile() = std::make_unique<H5Easy::File>(filename, H5Easy::File::Overwrite);
    LUPNT_CHECK(GetFile(), fmt::format("Failed to open log file: {}", filename.string()),
                "DataLogger");
  }

  // Log functions for different data types
  void DataLogger::Log(const std::string& name, double value) {
    std::lock_guard<std::mutex> lock(GetMutex());
    GetDoubleLogs()[name].push_back(value);
  }

  void DataLogger::Log(const std::string& name, Real value) {
    std::lock_guard<std::mutex> lock(GetMutex());
    GetDoubleLogs()[name].push_back(value.val());
  }

  void DataLogger::Log(const std::string& name, const std::string& value) {
    std::lock_guard<std::mutex> lock(GetMutex());
    GetStringLogs()[name].push_back(value);
  }

  void DataLogger::LogMatrix(const std::string& name, const Eigen::MatrixXd& value) {
    std::lock_guard<std::mutex> lock(GetMutex());
    GetMatrixLogs()[name].push_back(value);
  }

  void DataLogger::LogState(const std::string& parent, const State& state) {
    for (int i = 0; i < state.size(); i++) {
      auto name = state.GetNames()[i];
      auto unit = state.GetUnits()[i];
      std::replace(unit.begin(), unit.end(), '/', '_');
      auto value = state(i);
      Log(fmt::format("{}/{}_{}", parent, name, unit), value);
    }
  }

  std::string FindTimeKey(const std::string& name) {
    auto time_keys = std::vector<std::string>{"/time", "/t"};
    std::string target = "";

    // Parents
    std::string parent1 = name;
    size_t last_slash = parent1.rfind('/');
    if (last_slash != std::string::npos) parent1 = parent1.substr(0, last_slash);

    std::string parent2 = parent1;
    last_slash = parent2.rfind('/');
    if (last_slash != std::string::npos) parent2 = parent2.substr(0, last_slash);

    // Candidates
    bool found_time = false;
    std::vector<std::string> parents = {parent1, parent2};

    // Search
    for (const auto& [curr, data] : GetDoubleLogs()) {
      if (found_time) break;
      for (const auto& parent : parents) {
        for (const auto& time_key : time_keys) {
          std::string candidate = parent + time_key;
          if (curr == candidate) {
            target = candidate;
            found_time = true;
            break;
          }
        }
      }
    }
    return target;
  }

  // Flush all logs to HDF5 file
  void DataLogger::Flush() {
    if (!GetFile()) return;

    Logger::Info("Saving data", "DataLogger");
    std::lock_guard<std::mutex> lock(GetMutex());

    for (const auto& [name, data] : GetDoubleLogs()) {
      H5Easy::dump(*GetFile(), name, data);
    }
    for (const auto& [name, data] : GetStringLogs()) {
      H5Easy::dump(*GetFile(), name, data);
    }
    for (const auto& [name, data] : GetMatrixLogs()) {
      H5Easy::dump(*GetFile(), name, data);
    }

    GetFile()->flush();

    // Clear logs after writing
    GetDoubleLogs().clear();
    GetStringLogs().clear();
    GetMatrixLogs().clear();
  }

}  // namespace lupnt
