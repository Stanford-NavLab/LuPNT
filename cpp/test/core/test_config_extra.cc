#include <lupnt/core/config.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "../utils.cc"
#include "lupnt/devices/clock.h"

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("core.config_extra.read_enum") {
  // Integer form: enum_value indexes the enumerator list (ClockModel is 0-based
  // sequential: OCXO, USO, CSAC, ...), so index 1 -> USO, index 2 -> CSAC.
  YAML::Node n_int = YAML::Load("1");
  REQUIRE(ReadEnum<ClockModel>(n_int) == ClockModel::USO);

  YAML::Node n_int2 = YAML::Load("2");
  REQUIRE(ReadEnum<ClockModel>(n_int2) == ClockModel::CSAC);

  // String form goes through enum_cast by name.
  YAML::Node n_str = YAML::Load("RAFS");
  REQUIRE(ReadEnum<ClockModel>(n_str) == ClockModel::RAFS);

  // An unrecognized string fails both conversions and throws.
  YAML::Node n_bad = YAML::Load("__NOT_A_CLOCK__");
  REQUIRE_THROWS(ReadEnum<ClockModel>(n_bad));
}

TEST_CASE("core.config_extra.config_to_string") {
  Config config = YAML::Load("{alpha: 1, beta: two}");
  std::string s = ConfigToString(config);
  REQUIRE_THAT(s, ContainsSubstring("alpha"));
  REQUIRE_THAT(s, ContainsSubstring("beta"));
  REQUIRE_THAT(s, ContainsSubstring("two"));

  // Emitting then re-parsing preserves the values.
  Config reparsed = YAML::Load(s);
  REQUIRE(reparsed["alpha"].as<int>() == 1);
  REQUIRE(reparsed["beta"].as<std::string>() == "two");
}

TEST_CASE("core.config_extra.config_to_json_scalar_types") {
  Config config = YAML::Load(
      "{flag: true, i: 7, x: 3.5, name: hello, empty: null, seq: [10, 20], "
      "nested: {inner: 1}}");
  json j = ConfigToJson(config);

  REQUIRE(j["flag"].is_boolean());
  REQUIRE(j["flag"] == true);
  REQUIRE(j["i"].is_number_integer());
  REQUIRE(j["i"] == 7);
  REQUIRE(j["x"].is_number_float());
  REQUIRE(j["x"] == 3.5);
  REQUIRE(j["name"].is_string());
  REQUIRE(j["name"] == "hello");
  REQUIRE(j["empty"].is_null());
  REQUIRE(j["seq"].is_array());
  REQUIRE(j["seq"].size() == 2);
  REQUIRE(j["seq"][0] == 10);
  REQUIRE(j["nested"].is_object());
  REQUIRE(j["nested"]["inner"] == 1);

  // Round-trip back to a Config and confirm structure survives.
  Config back = JsonToConfig(j);
  REQUIRE(back["i"].as<int>() == 7);
  REQUIRE(back["nested"]["inner"].as<int>() == 1);
}

TEST_CASE("core.config_extra.load_config_node_overload") {
  // The Config-overload of LoadConfig clones and processes an in-memory node.
  Config in = YAML::Load("{a: 1, b: {c: 2}}");
  Config out = LoadConfig(in);
  REQUIRE(out["a"].as<int>() == 1);
  REQUIRE(out["b"]["c"].as<int>() == 2);
}

TEST_CASE("core.config_extra.load_config_auto_single_key") {
  auto dir = std::filesystem::temp_directory_path() / "lupnt_test_config_extra_auto";
  std::filesystem::create_directories(dir);
  auto path = dir / "wrapped.yaml";

  // A file whose root is a single-key map: LoadConfig with empty key
  // auto-extracts the inner map.
  Config config = YAML::Load("{wrapper: {x: 5, y: 6}}");
  SaveConfig(config, path.string());

  Config loaded = LoadConfig(path.string(), "", false);
  REQUIRE(loaded["x"].as<int>() == 5);
  REQUIRE(loaded["y"].as<int>() == 6);

  std::filesystem::remove_all(dir);
}

// NOTE: two further behaviors are intentionally left untested here:
//   * Local-key `inherit_from` via LoadConfig(path, key): specifying a key
//     makes LoadConfig reassign `node = root[key]`, which (yaml-cpp's
//     reference-assignment semantics) clobbers `root` so the sibling base key
//     is no longer visible to ProcessInheritance -- a pre-existing source
//     subtlety, not something a test should pin.
//   * Bare-filename lookup via config search dirs: FindConfigFile's first call
//     runs InitDefaultConfigSearchDirs -> GetBaseDir(), which requires the
//     LUPNT_PATH env var (not set under the test harness / CI), so it cannot be
//     exercised deterministically.
