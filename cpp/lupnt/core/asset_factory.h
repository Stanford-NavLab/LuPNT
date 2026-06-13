#pragma once
#include <functional>
#include <string>
#include <unordered_map>

#include "lupnt/core/config.h"
#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"

namespace lupnt {
  // Generalized Hybrid Registry/Factory Template
  // Usage: using MyFactory = AssetFactory<Base, Args...>;
  template <typename Base, typename... Args> class AssetFactory {
  public:
    using Creator = std::function<std::shared_ptr<Base>(Args...)>;
    static std::shared_ptr<Base> Create(const std::string& class_name, Args... args) {
      auto& reg = GetRegistry();
      auto it = reg.find(class_name);
      if (it != reg.end()) return it->second(std::forward<Args>(args)...);

      // Error
      std::stringstream ss;
      ss << "Unknown class: " << class_name << ". Available classes: ";
      for (const auto& entry : reg) ss << entry.first << " ";
      LUPNT_CHECK(false, ss.str(), "AssetFactory");
    }
    static void Register(const std::string& class_name, Creator creator) {
      Logger::Debug(fmt::format("Registering class: {}", class_name), "AssetFactory");
      GetRegistry()[class_name] = creator;
    }

    // Public accessor for registry (for debugging)
    static std::unordered_map<std::string, Creator>& GetRegistry();

  private:
    static std::unordered_map<std::string, Creator>& Registry() {
      static std::unordered_map<std::string, Creator> instance;
      return instance;
    }
  };

#define REGISTER_FACTORY_CLASS(BASE, DERIVED)                                               \
  namespace {                                                                               \
    struct DERIVED##Registrar {                                                             \
      DERIVED##Registrar() {                                                                \
        AssetFactory<BASE, Config&>::Register(#DERIVED,                                     \
                                              [](Config& config) -> std::shared_ptr<BASE> { \
                                                return std::make_shared<DERIVED>(config);   \
                                              });                                           \
      }                                                                                     \
    };                                                                                      \
    static DERIVED##Registrar global_##DERIVED##_registrar;                                 \
  }
}  // namespace lupnt
