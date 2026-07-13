#include <lupnt/core/config.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ---------------------------------------------------------------------------
// In-memory inheritance via LoadConfig(const Config&). This overload clones the
// node and runs ProcessInheritance against itself as the root, so local-key
// `inherit_from` references resolve without touching the filesystem or the
// config search dirs (which require LUPNT_PATH). Exercises the local-key
// branch and DeepMerge that test_config_extra.cc explicitly left uncovered.
// ---------------------------------------------------------------------------
TEST_CASE("core.config_more.local_inheritance_scalar_override") {
  Config in = YAML::Load(
      "base_profile: {x: 1, y: 2}\n"
      "active: {inherit_from: base_profile, y: 20, z: 30}\n");
  Config out = LoadConfig(in);

  // The base is preserved verbatim.
  REQUIRE(out["base_profile"]["x"].as<int>() == 1);
  REQUIRE(out["base_profile"]["y"].as<int>() == 2);

  // `active` inherits x from base, overrides y, and adds z. The inherit_from
  // marker itself is consumed.
  REQUIRE(out["active"]["x"].as<int>() == 1);
  REQUIRE(out["active"]["y"].as<int>() == 20);
  REQUIRE(out["active"]["z"].as<int>() == 30);
  REQUIRE_FALSE(out["active"]["inherit_from"]);
}

TEST_CASE("core.config_more.local_inheritance_nested_deep_merge") {
  // A nested map on the derived node must be deep-merged into the inherited
  // nested map (sub-keys combined, not wholesale replaced).
  Config in = YAML::Load(
      "base: {cfg: {a: 1, b: 2}}\n"
      "derived: {inherit_from: base, cfg: {b: 20, c: 3}}\n");
  Config out = LoadConfig(in);

  REQUIRE(out["derived"]["cfg"]["a"].as<int>() == 1);   // inherited untouched
  REQUIRE(out["derived"]["cfg"]["b"].as<int>() == 20);  // overridden
  REQUIRE(out["derived"]["cfg"]["c"].as<int>() == 3);   // added
}

TEST_CASE("core.config_more.local_inheritance_missing_key") {
  // A local inherit_from that names a key not present in the root falls back to
  // an empty map, into which the remaining keys are merged.
  Config in = YAML::Load("orphan: {inherit_from: __no_such_base__, a: 5}\n");
  Config out = LoadConfig(in);

  REQUIRE(out["orphan"]["a"].as<int>() == 5);
  REQUIRE_FALSE(out["orphan"]["inherit_from"]);
  // Nothing was inherited (the base key did not exist).
  REQUIRE_FALSE(out["orphan"]["x"]);
}

TEST_CASE("core.config_more.load_config_scalar_non_path_passthrough") {
  // LoadConfig(Config) on a scalar that is NOT a config path just clones/returns
  // it unchanged (no file lookup attempted).
  Config in = YAML::Load("just_a_word");
  Config out = LoadConfig(in);
  REQUIRE(out.IsScalar());
  REQUIRE(out.as<std::string>() == "just_a_word");
}

// ---------------------------------------------------------------------------
// File-based recursive loading. LoadRecursive resolves a scalar whose value is
// a *.yaml path into the loaded contents of that file. We use absolute paths so
// FindConfigFile short-circuits (is_absolute && exists) and never triggers the
// LUPNT_PATH-dependent default search-dir initialization.
// ---------------------------------------------------------------------------
TEST_CASE("core.config_more.load_recursive_absolute_subfile") {
  auto dir = std::filesystem::temp_directory_path() / "lupnt_test_config_more_rec";
  std::filesystem::create_directories(dir);
  auto sub = dir / "sub.yaml";
  auto main = dir / "main.yaml";

  Config sub_cfg = YAML::Load("{p: 10, q: 20}");
  SaveConfig(sub_cfg, sub.string());

  // Reference the sub-file by absolute path inside a scalar value.
  Config main_cfg;
  main_cfg["name"] = "root";
  main_cfg["sub"] = sub.string();
  SaveConfig(main_cfg, main.string());

  Config loaded = LoadConfig(main.string(), "", /*recursive=*/true);
  REQUIRE(loaded["name"].as<std::string>() == "root");
  // The scalar path was expanded into the loaded map.
  REQUIRE(loaded["sub"].IsMap());
  REQUIRE(loaded["sub"]["p"].as<int>() == 10);
  REQUIRE(loaded["sub"]["q"].as<int>() == 20);

  std::filesystem::remove_all(dir);
}

TEST_CASE("core.config_more.load_config_non_recursive_keeps_scalar_path") {
  // With recursive=false the *.yaml scalar is left as a plain string.
  auto dir = std::filesystem::temp_directory_path() / "lupnt_test_config_more_norec";
  std::filesystem::create_directories(dir);
  auto sub = dir / "sub.yaml";
  auto main = dir / "main.yaml";

  SaveConfig(YAML::Load("{p: 10}"), sub.string());
  Config main_cfg;
  // Two keys so the root is not a single-key map (which would auto-extract).
  main_cfg["name"] = "root";
  main_cfg["sub"] = sub.string();
  SaveConfig(main_cfg, main.string());

  Config loaded = LoadConfig(main.string(), "", /*recursive=*/false);
  REQUIRE(loaded["sub"].IsScalar());
  REQUIRE(loaded["sub"].as<std::string>() == sub.string());

  std::filesystem::remove_all(dir);
}

// ---------------------------------------------------------------------------
// ConfigToJson / JsonToConfig paths not covered by the map-of-scalars case in
// test_config_extra.cc: a top-level sequence of maps, and a bare top-level
// scalar going through the type-probing branch.
// ---------------------------------------------------------------------------
TEST_CASE("core.config_more.config_to_json_sequence_of_maps") {
  Config config = YAML::Load("[{id: 1, name: a}, {id: 2, name: b}]");
  json j = ConfigToJson(config);
  REQUIRE(j.is_array());
  REQUIRE(j.size() == 2);
  REQUIRE(j[0].is_object());
  REQUIRE(j[0]["id"] == 1);
  REQUIRE(j[0]["name"] == "a");
  REQUIRE(j[1]["id"] == 2);
  REQUIRE(j[1]["name"] == "b");

  // Round-trip back to a Config sequence.
  Config back = JsonToConfig(j);
  REQUIRE(back.IsSequence());
  REQUIRE(back.size() == 2);
  REQUIRE(back[1]["name"].as<std::string>() == "b");
}

TEST_CASE("core.config_more.config_to_json_top_level_scalars") {
  // Bool wins first, then int, then double, then string fallthrough.
  REQUIRE(ConfigToJson(YAML::Load("true")).is_boolean());
  REQUIRE(ConfigToJson(YAML::Load("true")) == true);
  REQUIRE(ConfigToJson(YAML::Load("42")).is_number_integer());
  REQUIRE(ConfigToJson(YAML::Load("42")) == 42);
  REQUIRE(ConfigToJson(YAML::Load("2.5")).is_number_float());
  REQUIRE(ConfigToJson(YAML::Load("2.5")) == 2.5);
  REQUIRE(ConfigToJson(YAML::Load("plain_string")).is_string());
  REQUIRE(ConfigToJson(YAML::Load("plain_string")) == "plain_string");
  REQUIRE(ConfigToJson(YAML::Load("~")).is_null());
}
