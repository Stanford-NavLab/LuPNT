#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "lupnt/core/config.h"
#include "py_yaml.h"  // NOLINT - needed for YAML::Node type caster

namespace py = pybind11;
using namespace lupnt;

// Forward declare ConfigWrapper so type caster can use it
class ConfigWrapper;

// Wrapper class for Config that provides attribute-style access
class ConfigWrapper {
public:
  Config node;

  ConfigWrapper() : node() {}
  explicit ConfigWrapper(const Config& n) : node(n) {}

  // Convert YAML node to Python, wrapping maps in ConfigWrapper
  static py::object node_to_python(const YAML::Node& n) {
    if (n.IsMap()) {
      return py::cast(ConfigWrapper(n));
    }
    return pybind11::detail::yaml_to_python(n);
  }

  // Get item with [] operator - raises KeyError if missing
  py::object get_item(const std::string& key) const {
    if (!node[key]) {
      throw py::key_error("Key '" + key + "' not found");
    }
    return node_to_python(node[key]);
  }

  // Set item with [] operator
  void set_item(const std::string& key, py::handle value) {
    node[key] = pybind11::detail::python_to_yaml(value);
  }

  // Get attribute with . operator - returns None if missing
  py::object get_attr(const std::string& key) const {
    if (!node[key]) {
      return py::none();
    }
    return node_to_python(node[key]);
  }

  // Contains check
  bool contains(const std::string& key) const { return node[key].IsDefined(); }

  // Size
  size_t size() const { return node.size(); }

  // Iterate over keys
  py::list keys() const {
    py::list result;
    for (auto it = node.begin(); it != node.end(); ++it) {
      result.append(it->first.as<std::string>());
    }
    return result;
  }

  // Iterate over values
  py::list values() const {
    py::list result;
    for (auto it = node.begin(); it != node.end(); ++it) {
      result.append(node_to_python(it->second));
    }
    return result;
  }

  // Iterate over items
  py::list items() const {
    py::list result;
    for (auto it = node.begin(); it != node.end(); ++it) {
      py::tuple item = py::make_tuple(it->first.as<std::string>(), node_to_python(it->second));
      result.append(item);
    }
    return result;
  }

  // Get with default
  py::object get(const std::string& key, py::object default_value = py::none()) const {
    if (!node[key]) {
      return default_value;
    }
    return node_to_python(node[key]);
  }

  // Convert to string
  std::string to_string() const { return ConfigToString(node); }

  // Load from dict
  void update(py::dict dict) { node = pybind11::detail::python_to_yaml(dict); }

  // To dict
  py::dict to_dict() const { return pybind11::detail::yaml_to_python(node).cast<py::dict>(); }

  // Get the underlying YAML::Node
  const YAML::Node& get_node() const { return node; }

  // Copy support
  ConfigWrapper copy() const { return ConfigWrapper(YAML::Clone(node)); }

  // Pickle support - return state as dict
  py::dict get_state() const { return to_dict(); }

  // Pickle support - restore from dict
  void set_state(py::dict state) { node = pybind11::detail::python_to_yaml(state); }
};

void InitConfig(py::module& m) {
  // Config wrapper class
  auto config_class
      = py::class_<ConfigWrapper>(m, "Config")
            .def(py::init<>())
            .def(py::init<const Config&>())
            .def(py::init([](py::dict dict) {
              Config config = pybind11::detail::python_to_yaml(dict);
              return ConfigWrapper(config);
            }))
            .def("__getitem__", &ConfigWrapper::get_item)
            .def("__setitem__", &ConfigWrapper::set_item)
            .def("__getattr__", &ConfigWrapper::get_attr)
            .def("__contains__", &ConfigWrapper::contains)
            .def("__len__", &ConfigWrapper::size)
            .def("keys", &ConfigWrapper::keys)
            .def("values", &ConfigWrapper::values)
            .def("items", &ConfigWrapper::items)
            .def("get", &ConfigWrapper::get, py::arg("key"), py::arg("default") = py::none())
            .def("to_string", &ConfigWrapper::to_string)
            .def("update", &ConfigWrapper::update)
            .def("to_dict", &ConfigWrapper::to_dict)
            .def("__repr__", &ConfigWrapper::to_string)
            .def("__str__", &ConfigWrapper::to_string)
            // Copy support
            .def("copy", &ConfigWrapper::copy, "Create a deep copy of the config")
            .def("__copy__", [](const ConfigWrapper& self) { return self.copy(); })
            .def(
                "__deepcopy__", [](const ConfigWrapper& self, py::dict) { return self.copy(); },
                py::arg("memo"))
            // Pickle support
            .def(py::pickle(
                [](const ConfigWrapper& self) {  // __getstate__
                  return self.get_state();
                },
                [](py::dict state) {  // __setstate__
                  ConfigWrapper config;
                  config.set_state(state);
                  return config;
                }));

  // Config management functions
  m.def("add_config_search_dir", &AddConfigSearchDir, py::arg("dir"));
  m.def("init_default_config_search_dirs", &InitDefaultConfigSearchDirs);

  // load_config with path (string)
  m.def(
      "load_config",
      [](const std::string& path, const std::string& key, bool recursive) {
        Config config = LoadConfig(path, key, recursive);
        return ConfigWrapper(config);
      },
      py::arg("path"), py::arg("key") = "", py::arg("recursive") = true);

  // load_config with Config (accepts dict or YAML::Node via type caster)
  m.def(
      "load_config",
      [](const Config& config) {
        Config result = LoadConfig(config);
        return ConfigWrapper(result);
      },
      py::arg("config"));

  // load_config with ConfigWrapper
  m.def(
      "load_config",
      [](const ConfigWrapper& config_wrapper) {
        Config result = LoadConfig(config_wrapper.node);
        return ConfigWrapper(result);
      },
      py::arg("config"));

  m.def(
      "save_config",
      [](const ConfigWrapper& config_wrapper, const std::string& path) {
        SaveConfig(config_wrapper.node, path);
      },
      py::arg("config"), py::arg("path"));

  m.def(
      "config_to_string",
      [](const ConfigWrapper& config_wrapper) { return ConfigToString(config_wrapper.node); },
      py::arg("config"));
}
