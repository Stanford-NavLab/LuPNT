#include "lupnt/core/config.h"

#include <algorithm>
#include <fstream>

#include "lupnt/core/error.h"
#include "lupnt/core/file.h"

namespace lupnt {

  static std::vector<std::filesystem::path> config_search_dirs_;

  void AddConfigSearchDir(const std::string& dir) {
    std::filesystem::path p(dir);
    if (std::filesystem::exists(p) && std::filesystem::is_directory(p)
        && std::find(config_search_dirs_.begin(), config_search_dirs_.end(), p)
               == config_search_dirs_.end()) {
      config_search_dirs_.push_back(p);
    }
  }

  void InitDefaultConfigSearchDirs() { AddConfigSearchDir((GetBaseDir() / "configs").string()); }

  // Parse "path::key" into (path, key)
  static std::pair<std::string, std::string> ParsePathAndKey(const std::string& path_with_key) {
    size_t pos = path_with_key.find("::");
    return pos != std::string::npos
               ? std::make_pair(path_with_key.substr(0, pos), path_with_key.substr(pos + 2))
               : std::make_pair(path_with_key, "");
  }

  static bool IsConfigPath(const std::string& str) {
    return str.find(".yaml") != std::string::npos || str.find(".yml") != std::string::npos;
  }

  static std::filesystem::path FindConfigFile(const std::string& path);
  static void LoadRecursive(YAML::Node& node);
  static YAML::Node ProcessInheritance(const YAML::Node& node, const YAML::Node& root);

  static YAML::Node ProcessScalarValue(const YAML::Node& value) {
    if (!value.IsScalar()) return value;
    std::string str = value.as<std::string>();
    if (!IsConfigPath(str)) return value;

    try {
      auto [path, key] = ParsePathAndKey(str);
      auto abs_path = FindConfigFile(path);
      if (std::filesystem::exists(abs_path)) {
        return LoadConfig(abs_path.string(), key, true);
      }
    } catch (...) {
    }
    return value;
  }

  static void DeepMerge(YAML::Node& base, const std::string& key, const YAML::Node& value) {
    if (!(value.IsMap() && base[key] && base[key].IsMap())) {
      base[key] = ProcessScalarValue(value);
      return;
    }

    for (auto it = value.begin(); it != value.end(); ++it) {
      std::string sub_key = it->first.as<std::string>();
      if (base[key][sub_key] && base[key][sub_key].IsMap() && it->second.IsMap()) {
        for (auto sub_it = it->second.begin(); sub_it != it->second.end(); ++sub_it) {
          base[key][sub_key][sub_it->first.as<std::string>()] = ProcessScalarValue(sub_it->second);
        }
      } else {
        base[key][sub_key] = ProcessScalarValue(it->second);
      }
    }
  }

  static YAML::Node ProcessInheritance(const YAML::Node& node, const YAML::Node& root) {
    if (!node.IsMap()) return node;

    YAML::Node result;

    // First, handle inheritance at this level if present
    if (node["inherit_from"]) {
      std::string inherit_path = node["inherit_from"].as<std::string>();
      try {
        // Check if this is a simple key reference (no file path)
        if (!IsConfigPath(inherit_path) && inherit_path.find("::") == std::string::npos) {
          // Try to find the key in the root node of the current file
          if (root[inherit_path]) {
            result = YAML::Clone(root[inherit_path]);
            Logger::Debug(fmt::format("Inherited config from local key: {}", inherit_path),
                          "Config");
          } else {
            Logger::Warn(fmt::format("Local key '{}' not found in config", inherit_path), "Config");
            result = YAML::Node(YAML::NodeType::Map);
          }
        } else {
          // Load from external file
          auto [file, key] = ParsePathAndKey(inherit_path);
          result = LoadConfig(file, key, true);
        }

        for (auto it = node.begin(); it != node.end(); ++it) {
          if (it->first.as<std::string>() != "inherit_from") {
            DeepMerge(result, it->first.as<std::string>(), it->second);
          }
        }
      } catch (const std::exception& e) {
        Logger::Warn(fmt::format("Failed to inherit from {}: {}", inherit_path, e.what()),
                     "Config");
        result = YAML::Clone(node);
      }
    } else {
      result = YAML::Clone(node);
    }

    // Recursively process inheritance for all nested maps
    for (auto it = result.begin(); it != result.end(); ++it) {
      if (it->second.IsMap()) {
        it->second = ProcessInheritance(it->second, root);
      }
    }

    return result;
  }

  Config LoadConfig(const std::string& path_with_key, const std::string& key, bool recursive) {
    auto [parsed_path, parsed_key] = ParsePathAndKey(path_with_key);
    std::string actual_key = key.empty() ? parsed_key : key;

    auto abs_path = FindConfigFile(parsed_path);
    Logger::Info(fmt::format("Loading config file: {}", abs_path.string()), "Config");
    YAML::Node root = YAML::LoadFile(abs_path.string());

    // Add parent directory to search dirs
    auto parent = abs_path.parent_path();
    if (std::find(config_search_dirs_.begin(), config_search_dirs_.end(), parent)
        == config_search_dirs_.end()) {
      config_search_dirs_.push_back(parent);
    }

    // Extract key or auto-extract single-key map
    YAML::Node node = root;
    if (!actual_key.empty()) {
      node = root[actual_key];
    } else if (root.IsMap() && root.size() == 1) {
      node = root.begin()->second;
    }

    if (recursive) LoadRecursive(node);
    return ProcessInheritance(node, root);
  }

  Config LoadConfig(const Config& config) {
    // If config is just a string that looks like a file path, load that file
    if (config.IsScalar()) {
      std::string str = config.as<std::string>();
      if (IsConfigPath(str)) {
        return LoadConfig(str, "", true);
      }
    }
    YAML::Node cloned = YAML::Clone(config);
    return ProcessInheritance(cloned, cloned);
  }

  std::string ConfigToString(const Config& config) {
    YAML::Emitter out;
    out << config;
    return out.c_str();
  }

  void SaveConfig(const Config& config, const std::string& path) { std::ofstream(path) << config; }

  static bool PathEndsWith(const std::filesystem::path& rel_path,
                           const std::filesystem::path& req_path) {
    std::vector<std::filesystem::path> rel, req;
    for (const auto& c : rel_path) rel.push_back(c);
    for (const auto& c : req_path) req.push_back(c);

    if (req.size() > rel.size()) return false;
    for (size_t i = 0; i < req.size(); ++i) {
      if (rel[rel.size() - req.size() + i] != req[i]) return false;
    }
    return true;
  }

  static std::filesystem::path FindConfigFile(const std::string& path) {
    std::filesystem::path p(path);
    if (p.is_absolute() && std::filesystem::exists(p)) return p;

    static bool initialized = false;
    if (!initialized) {
      InitDefaultConfigSearchDirs();
      initialized = true;
    }

    Logger::Debug(fmt::format("Searching for config file: {}", path), "Config");

    for (auto it = config_search_dirs_.rbegin(); it != config_search_dirs_.rend(); ++it) {
      if (std::filesystem::exists(*it / p)) return *it / p;

      try {
        for (const auto& entry : std::filesystem::recursive_directory_iterator(*it)) {
          if (entry.is_regular_file() && PathEndsWith(entry.path().lexically_relative(*it), p)) {
            Logger::Debug(fmt::format("Found: {}", entry.path().string()), "Config");
            return entry.path();
          }
        }
      } catch (const std::filesystem::filesystem_error&) {
      }
    }

    if (std::filesystem::exists(p)) return p;
    LUPNT_CHECK(false, "Config file not found: " + path, "Config");
  }

  static void LoadRecursive(YAML::Node& node) {
    if (!node.IsMap()) return;

    for (auto it = node.begin(); it != node.end(); ++it) {
      if (it->second.IsScalar()) {
        std::string str = it->second.as<std::string>();
        if (IsConfigPath(str)) {
          try {
            auto [path, key] = ParsePathAndKey(str);
            auto abs_path = FindConfigFile(path);
            if (std::filesystem::exists(abs_path)) {
              it->second = LoadConfig(abs_path.string(), key, true);
            }
          } catch (...) {
          }
        }
      } else if (it->second.IsMap()) {
        LoadRecursive(it->second);
      }
    }
  }

  json ConfigToJson(const Config& config) {
    if (config.IsNull()) {
      return nullptr;
    } else if (config.IsScalar()) {
      if (config.Tag() == "!") {
        return config.as<std::string>();
      }
      // Try to parse as different types
      try {
        return config.as<bool>();
      } catch (...) {
      }
      try {
        return config.as<int>();
      } catch (...) {
      }
      try {
        return config.as<double>();
      } catch (...) {
      }
      return config.as<std::string>();
    } else if (config.IsSequence()) {
      json arr = json::array();
      for (const auto& item : config) {
        arr.push_back(ConfigToJson(item));
      }
      return arr;
    } else if (config.IsMap()) {
      json obj = json::object();
      for (const auto& pair : config) {
        obj[pair.first.as<std::string>()] = ConfigToJson(pair.second);
      }
      return obj;
    }
    return nullptr;
  }

  Config JsonToConfig(const json& json) { return YAML::Load(json.dump()); }
}  // namespace lupnt
