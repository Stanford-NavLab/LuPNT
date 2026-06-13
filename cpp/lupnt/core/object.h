#pragma once

#include <atomic>
#include <mutex>
#include <string>
#include <typeinfo>
#include <unordered_map>

namespace lupnt {

  template <typename Derived> class Object {
  private:
    static std::unordered_map<std::string, std::atomic<int>>& GetInstanceMap() {
      static std::unordered_map<std::string, std::atomic<int>> instance_counts;
      return instance_counts;
    }

    static std::mutex& GetMutex() {
      static std::mutex mutex;
      return mutex;
    }

    std::string id_;

    // Extract class name **without** namespace
    static std::string ExtractClassName(const char* pretty_func) {
      std::string class_name(pretty_func);

#if defined(__clang__) || defined(__GNUC__)
      size_t start = class_name.find("Derived = ") + 10;
      size_t end = class_name.find("]", start);
#elif defined(_MSC_VER)
      size_t start = class_name.find("Object<") + 7;
      size_t end = class_name.find(">", start);
#endif

      if (start != std::string::npos && end != std::string::npos) {
        class_name = class_name.substr(start, end - start);
      } else {
        class_name = "Unknown";
      }

      // Remove namespace prefix (lupnt:: or any other)
      size_t ns_pos = class_name.rfind("::");
      if (ns_pos != std::string::npos) {
        class_name = class_name.substr(ns_pos + 2);
      }

      return class_name;
    }

  protected:
    void InitializeId() {
      std::string class_name = ExtractClassName(__PRETTY_FUNCTION__);

      std::lock_guard<std::mutex> lock(GetMutex());
      int instance_number = GetInstanceMap()[class_name]++;
      id_ = class_name + "_" + std::to_string(instance_number);
    }

  public:
    Object() { InitializeId(); }
    virtual ~Object() = default;

    std::string GetId() const { return id_; }
    void SetId(const std::string& id) { id_ = id; }

    static int GetInstanceCount(const std::string& class_name) {
      return GetInstanceMap()[class_name];
    }
  };

}  // namespace lupnt
