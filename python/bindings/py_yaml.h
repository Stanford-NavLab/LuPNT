#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <yaml-cpp/yaml.h>

namespace py = pybind11;

namespace pybind11 {
  namespace detail {

    // Convert Python object to YAML::Node (used by multiple type casters)
    inline YAML::Node python_to_yaml(handle obj) {
      if (obj.is_none()) {
        return YAML::Node();
      } else if (py::isinstance<py::bool_>(obj)) {
        return YAML::Node(obj.cast<bool>());
      } else if (py::isinstance<py::int_>(obj)) {
        return YAML::Node(obj.cast<int64_t>());
      } else if (py::isinstance<py::float_>(obj)) {
        return YAML::Node(obj.cast<double>());
      } else if (py::isinstance<py::str>(obj)) {
        return YAML::Node(obj.cast<std::string>());
      } else if (py::isinstance<py::list>(obj) || py::isinstance<py::tuple>(obj)) {
        YAML::Node node(YAML::NodeType::Sequence);
        for (auto item : obj) {
          node.push_back(python_to_yaml(item));
        }
        return node;
      } else if (py::isinstance<py::dict>(obj)) {
        YAML::Node node(YAML::NodeType::Map);
        for (auto item : py::cast<py::dict>(obj)) {
          std::string key = py::str(item.first);
          node[key] = python_to_yaml(item.second);
        }
        return node;
      } else if (py::hasattr(obj, "__fspath__")) {
        // Handle pathlib.Path and other path-like objects
        return YAML::Node(py::str(obj).cast<std::string>());
      } else if (py::hasattr(obj, "__array__") || py::hasattr(obj, "__iter__")) {
        // Handle numpy arrays and other iterables (but not strings, already handled above)
        try {
          YAML::Node node(YAML::NodeType::Sequence);
          for (auto item : obj) {
            node.push_back(python_to_yaml(item));
          }
          return node;
        } catch (...) {
          throw std::runtime_error("Unsupported Python type for YAML conversion");
        }
      } else {
        throw std::runtime_error("Unsupported Python type for YAML conversion");
      }
    }

    // Convert YAML::Node to Python object (used by multiple type casters)
    inline py::object yaml_to_python(const YAML::Node& node) {
      switch (node.Type()) {
        case YAML::NodeType::Null: return py::none();

        case YAML::NodeType::Scalar: {
          // Try to infer type
          try {
            return py::cast(node.as<bool>());
          } catch (...) {
          }
          try {
            return py::cast(node.as<int64_t>());
          } catch (...) {
          }
          try {
            return py::cast(node.as<double>());
          } catch (...) {
          }
          return py::cast(node.as<std::string>());
        }

        case YAML::NodeType::Sequence: {
          py::list result;
          for (const auto& item : node) {
            result.append(yaml_to_python(item));
          }
          return result;
        }

        case YAML::NodeType::Map: {
          py::dict result;
          for (const auto& item : node) {
            std::string key = item.first.as<std::string>();
            result[py::cast(key)] = yaml_to_python(item.second);
          }
          return result;
        }

        default: return py::none();
      }
    }

    // Type caster for YAML::Node <-> Python dict/list/scalar
    template <> struct type_caster<YAML::Node> {
    public:
      PYBIND11_TYPE_CASTER(YAML::Node, _("dict"));

      // Python -> C++ conversion
      bool load(handle src, bool) {
        if (src.is_none()) {
          value = YAML::Node();
          return true;
        }

        // Check if it's a ConfigWrapper by checking for get_node method
        // ConfigWrapper has a public get_node() method that returns the YAML::Node
        if (hasattr(src, "keys") && hasattr(src, "to_dict")) {
          try {
            // It quacks like a ConfigWrapper - try to convert via to_dict()
            auto dict = src.attr("to_dict")();
            value = python_to_yaml(dict);
            return true;
          } catch (...) {
            // Fall through to regular conversion
          }
        }

        try {
          value = python_to_yaml(src);
          return true;
        } catch (...) {
          return false;
        }
      }

      // C++ -> Python conversion
      static handle cast(const YAML::Node& src, return_value_policy, handle) {
        return yaml_to_python(src).release();
      }
    };

  }  // namespace detail
}  // namespace pybind11

// No need for a separate Config type caster since Config is just an alias for YAML::Node
