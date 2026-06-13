#pragma once

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace pecsim::examples {

  struct ExampleConfig {
    int year = 2025;
    int month = 3;
    int day = 1;
    int hour = 12;
    int minute = 0;
    double second = 0.0;
    double kp = 3.0;
    double rz12 = 50.0;
    bool run_cpp = true;
    bool run_fortran = false;
    bool only_xz = true;
  };

  inline std::string trim_copy(std::string value) {
    auto not_space = [](unsigned char ch) { return !std::isspace(ch); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
  }

  inline bool parse_bool(const std::string& value) {
    std::string lowered = value;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char ch) { return std::tolower(ch); });

    if (lowered == "1" || lowered == "true" || lowered == "yes" || lowered == "on") {
      return true;
    }
    if (lowered == "0" || lowered == "false" || lowered == "no" || lowered == "off") {
      return false;
    }
    throw std::runtime_error("Invalid boolean value: " + value);
  }

  inline ExampleConfig load_example_config(const std::filesystem::path& config_path) {
    std::ifstream input(config_path);
    if (!input.is_open()) {
      throw std::runtime_error("Could not open config file: " + config_path.string());
    }

    ExampleConfig config;
    std::string line;
    int line_number = 0;

    while (std::getline(input, line)) {
      ++line_number;
      const std::string trimmed = trim_copy(line);
      if (trimmed.empty() || trimmed[0] == '#') {
        continue;
      }

      const size_t eq = trimmed.find('=');
      if (eq == std::string::npos) {
        throw std::runtime_error("Invalid config entry at line " + std::to_string(line_number)
                                 + ": " + trimmed);
      }

      const std::string key = trim_copy(trimmed.substr(0, eq));
      const std::string value = trim_copy(trimmed.substr(eq + 1));

      if (key == "year") {
        config.year = std::stoi(value);
      } else if (key == "month") {
        config.month = std::stoi(value);
      } else if (key == "day") {
        config.day = std::stoi(value);
      } else if (key == "hour") {
        config.hour = std::stoi(value);
      } else if (key == "minute") {
        config.minute = std::stoi(value);
      } else if (key == "second") {
        config.second = std::stod(value);
      } else if (key == "kp") {
        config.kp = std::stod(value);
      } else if (key == "rz12") {
        config.rz12 = std::stod(value);
      } else if (key == "run_cpp") {
        config.run_cpp = parse_bool(value);
      } else if (key == "run_fortran") {
        config.run_fortran = parse_bool(value);
      } else if (key == "only_xz") {
        config.only_xz = parse_bool(value);
      } else {
        throw std::runtime_error("Unknown config key at line " + std::to_string(line_number) + ": "
                                 + key);
      }
    }

    return config;
  }

}  // namespace pecsim::examples
