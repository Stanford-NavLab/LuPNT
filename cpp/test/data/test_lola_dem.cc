#include <lupnt/interfaces/lola_dem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("data.lola_dem.sites") {
  const auto& sites = GetLolaSites();
  REQUIRE(sites.size() == 6);

  // Query near the Connecting Ridge center selects Site01.
  REQUIRE(SelectLolaSite(-89.45, 222.8).id == "Site01");
  // Query near the Malapert massif selects Site23.
  REQUIRE(SelectLolaSite(-86.0, 3.0).id == "Site23");
  // Query near the Leibnitz beta plateau selects Site20.
  REQUIRE(SelectLolaSite(-85.3, 33.0).id == "Site20");

  REQUIRE(LolaDemUrl("Site01")
          == "https://pgda.gsfc.nasa.gov/data/LOLA_5mpp/Site01/Site01_final_adj_5mpp_surf.tif");
}

TEST_CASE("data.lola_dem.bilinear") {
  // 3x3 regular grid: 10 m spacing, x increases with column, y decreases with row.
  const int ny = 3, nx = 3;
  MatXd x(ny, nx), y(ny, nx), elev(ny, nx);
  for (int r = 0; r < ny; ++r) {
    for (int c = 0; c < nx; ++c) {
      x(r, c) = 100.0 + c * 10.0;
      y(r, c) = 200.0 - r * 10.0;
      elev(r, c) = c * 1.0 + r * 2.0;  // planar field
    }
  }
  LolaSite site{"SiteTest", "test", -89.0, 0.0};
  LunarDem dem(x, y, elev, site);

  // Exact grid points.
  REQUIRE_THAT(dem.GetElevation(100.0, 200.0), WithinAbs(0.0, 1e-9));
  REQUIRE_THAT(dem.GetElevation(120.0, 200.0), WithinAbs(2.0, 1e-9));
  REQUIRE_THAT(dem.GetElevation(100.0, 180.0), WithinAbs(4.0, 1e-9));
  // Bilinear midpoint (col 0.5 -> +0.5, row 0.5 -> +1.0).
  REQUIRE_THAT(dem.GetElevation(105.0, 195.0), WithinAbs(1.5, 1e-9));
  // Out-of-range queries clamp to the nearest edge.
  REQUIRE_THAT(dem.GetElevation(0.0, 1000.0), WithinAbs(0.0, 1e-9));

  REQUIRE(dem.rows() == 3);
  REQUIRE(dem.cols() == 3);
  REQUIRE_THAT(dem.center_x(), WithinAbs(110.0, 1e-9));
  REQUIRE_THAT(dem.center_y(), WithinAbs(190.0, 1e-9));
}
