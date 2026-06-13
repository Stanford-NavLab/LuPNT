#include <gdal_priv.h>
#include <lupnt/data/dem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("data.dem") {
  GDALAllRegister();
  std::filesystem::path path = std::filesystem::temp_directory_path() / "lupnt_dem_test.tif";
  GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
  REQUIRE(driver != nullptr);

  GDALDataset* dataset = driver->Create(path.string().c_str(), 4, 3, 1, GDT_Float32, nullptr);
  REQUIRE(dataset != nullptr);
  double geotransform[6] = {100.0, 10.0, 0.0, 200.0, 0.0, -10.0};
  REQUIRE(dataset->SetGeoTransform(geotransform) == CE_None);

  std::vector<float> values{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  REQUIRE(dataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, 4, 3, values.data(), 4, 3,
                                              GDT_Float32, 0, 0)
          == CE_None);
  GDALClose(dataset);

  auto [x, y, data] = LoadTiff(path, {110.0, 130.0}, {170.0, 190.0}, 10.0);
  REQUIRE(x.rows() == 2);
  REQUIRE(x.cols() == 2);
  REQUIRE(y.rows() == 2);
  REQUIRE(data.rows() == 2);
  REQUIRE(data.cols() == 2);
  REQUIRE_THAT(x(0, 0), WithinAbs(115.0, epsilon));
  REQUIRE_THAT(y(0, 0), WithinAbs(185.0, epsilon));
  REQUIRE_THAT(data(0, 0), WithinAbs(5.0, epsilon));
  REQUIRE_THAT(data(1, 1), WithinAbs(10.0, epsilon));
}
