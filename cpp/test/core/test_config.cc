#include <lupnt/core/config.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.config") {
  SECTION("YAML config converts to JSON and back") {
    Config config = YAML::Load("{name: demo, count: 3, enabled: true, values: [1, 2]}");

    json encoded = ConfigToJson(config);
    Config decoded = JsonToConfig(encoded);

    REQUIRE(encoded["name"] == "demo");
    REQUIRE(encoded["count"] == 3);
    REQUIRE(encoded["enabled"] == true);
    REQUIRE(decoded["values"][1].as<int>() == 2);
  }

  SECTION("SaveConfig and LoadConfig round-trip a selected key") {
    auto dir = std::filesystem::temp_directory_path() / "lupnt_test_config";
    std::filesystem::create_directories(dir);
    auto path = dir / "config.yaml";

    Config config = YAML::Load("{root: {answer: 42, label: ok}}");
    SaveConfig(config, path.string());

    Config loaded = LoadConfig(path.string(), "root", false);

    REQUIRE(loaded["answer"].as<int>() == 42);
    REQUIRE(loaded["label"].as<std::string>() == "ok");

    std::filesystem::remove_all(dir);
  }
}
