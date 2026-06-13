#pragma once

#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <filesystem>
#include <highfive/H5Easy.hpp>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/states/state.h"

#define LOG(id, ...) lupnt::DataLogger::LogHelper(id, #__VA_ARGS__, __VA_ARGS__)

namespace lupnt {

  class DataLogger {
  public:
    using MatrixMap = std::unordered_map<std::string, std::vector<Eigen::MatrixXd>>;
    using DoubleMap = std::unordered_map<std::string, std::vector<double>>;
    using StringMap = std::unordered_map<std::string, std::vector<std::string>>;

  public:
    static void SetOutputFile(const std::filesystem::path& filename);

    static void Log(const std::string& name, double value);
    static void Log(const std::string& name, Real value);
    static void Log(const std::string& name, const std::string& value);

    static void LogState(const std::string& parent, const State& state);

    template <typename Derived>
    static void Log(const std::string& name, const Eigen::DenseBase<Derived>& value) {
      // Forward to the implementation in the .cc file
      LogMatrix(name, value.template cast<double>());
    }

    static void Flush();

    template <typename... Args>
    static void LogHelper(const std::string& id, const std::string& names, const Args&... args) {
      std::istringstream name_stream(names);
      std::vector<std::string> name_list;
      std::string name;

      // Efficiently split the comma-separated names
      while (std::getline(name_stream, name, ',')) {
        // Trim leading/trailing spaces
        name.erase(0, name.find_first_not_of(" \t"));
        name.erase(name.find_last_not_of(" \t") + 1);
        // Remove final _ if present
        if (!name.empty() && name.back() == '_') name.pop_back();
        name_list.push_back(name);
      }

      // Ensure number of names matches number of arguments
      constexpr size_t num_args = sizeof...(args);
      if (name_list.size() != num_args) {
        spdlog::error("Mismatch: Expected {} names but received {} arguments.", name_list.size(),
                      num_args);
        return;
      }

      // Log each variable with its corresponding name, prefixed by the ID
      size_t index = 0;
      ((Log(id + "/" + name_list[index++], args)), ...);
    }

    virtual void Log(Real time) = 0;

  private:
    static void LogMatrix(const std::string& name, const Eigen::MatrixXd& value);
  };

}  // namespace lupnt
