#include <lupnt/core/user_file_path.h>
#include <lupnt/datasets/crater_data.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdlib>
#include <filesystem>
#include <iostream>

using namespace lupnt;
TEST_CASE("CraterDataLoader Test", "[CraterDataLoader]") {
  std::string filename = "lunar_crater_database_robbins_2018.csv";
  std::filesystem::path filepath = GetFilePath(filename);

  REQUIRE(std::filesystem::exists(filepath));

  SECTION("Load Craters") {
    try {
      std::vector<Crater> craters = CraterDataLoader::LoadCraters(filename);
      REQUIRE(!craters.empty());

      // Check some properties of the loaded data (this is an example, adjust
      // based on your actual data)
      REQUIRE(craters[0].lat >= -90.0);
      REQUIRE(craters[0].lat <= 90.0);
      REQUIRE(craters[0].lon >= -180.0);
      REQUIRE(craters[0].lon <= 180.0);
      REQUIRE(craters[0].diam_major >= 0.0);
      REQUIRE(craters[0].diam_minor >= 0.0);
    } catch (const std::exception& e) {
      FAIL("Exception occurred: " << e.what());
    }
  }

  SECTION("Extract Robbins Dataset") {
    try {
      std::vector<Crater> craters = CraterDataLoader::LoadCraters(filename);
      REQUIRE(!craters.empty());

      Eigen::VectorXd lat, lon, major, minor, psi;
      std::vector<std::string> crater_id;
      CraterDataLoader::ExtractRobbinsDataset(craters, lat, lon, major, minor,
                                              psi, crater_id);

      REQUIRE(lat.size() == craters.size());
      REQUIRE(lon.size() == craters.size());
      REQUIRE(major.size() == craters.size());
      REQUIRE(minor.size() == craters.size());
      REQUIRE(psi.size() == craters.size());
      REQUIRE(crater_id.size() == craters.size());

      // Check some properties of the extracted data (this is an example, adjust
      // based on your actual data)
      REQUIRE(lat[0] >= -M_PI / 2);
      REQUIRE(lat[0] <= M_PI / 2);
      REQUIRE(lon[0] >= -M_PI);
      REQUIRE(lon[0] <= M_PI);
    } catch (const std::exception& e) {
      FAIL("Exception occurred: " << e.what());
    }
  }
}