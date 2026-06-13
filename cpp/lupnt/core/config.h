#pragma once

#include <yaml-cpp/yaml.h>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>
#include <string>

#include "lupnt/core/definitions.h"

namespace lupnt {

  // Config loading and management functions
  void AddConfigSearchDir(const std::string& dir);
  void InitDefaultConfigSearchDirs();

  Config LoadConfig(const std::string& yaml_path, const std::string& key = "",
                    bool recursive = true);
  Config LoadConfig(const Config& config);

  void SaveConfig(const Config& config, const std::string& path);

  std::string ConfigToString(const Config& config);

  json ConfigToJson(const Config& config);
  Config JsonToConfig(const json& json);

  // Reads an enum value from a YAML node, accepting either an integer or a string.
  // Throws if both conversions fail.
  template <typename T> T ReadEnum(const YAML::Node& node) {
    try {
      // Try to read as int (numeric enum value)
      return enum_value<T>(node.as<int>());
    } catch (const YAML::BadConversion&) {
      // Otherwise try string
      auto result = enum_cast<T>(node.as<std::string>());
      if (result.has_value())
        return result.value();
      else
        throw YAML::BadConversion(node.Mark());
    }
  }

}  // namespace lupnt
