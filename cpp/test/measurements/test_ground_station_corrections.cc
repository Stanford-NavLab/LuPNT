#include <lupnt/core/constants.h>
#include <lupnt/measurements/ground_station_corrections.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using namespace lupnt;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("ground_station_corrections.troposphere") {
  // Saastamoinen zenith total delay at sea level is ~2.3-2.6 m.
  double zenith = TroposphereDelaySaastamoinen(PI / 2, 35.0 * RAD, 0.0, 1013.25, 288.15, 50.0);
  REQUIRE(zenith > 2.2);
  REQUIRE(zenith < 2.8);

  // The obliquity mapping is 1/sin(el): the slant delay at 10 deg equals zenith / sin(10 deg).
  double el10 = TroposphereDelaySaastamoinen(10.0 * RAD, 35.0 * RAD, 0.0, 1013.25, 288.15, 50.0);
  REQUIRE_THAT(el10, WithinRel(zenith / std::sin(10.0 * RAD), 1e-9));
  REQUIRE(el10 > zenith);

  // A drier atmosphere reduces the wet delay, so the total shrinks.
  double dry = TroposphereDelaySaastamoinen(PI / 2, 35.0 * RAD, 0.0, 1013.25, 288.15, 0.0);
  REQUIRE(dry < zenith);

  // Standard-atmosphere helpers decrease with height.
  REQUIRE(StandardAtmospherePressureHPa(0.0) > StandardAtmospherePressureHPa(1000.0));
  REQUIRE(StandardAtmosphereTemperatureK(0.0) > StandardAtmosphereTemperatureK(1000.0));
  REQUIRE_THAT(StandardAtmospherePressureHPa(0.0), WithinAbs(1013.25, 1e-6));
}

TEST_CASE("ground_station_corrections.ionosphere") {
  // Delay = 40.3 * STEC / f^2. At 10 TECU, X-band (8.4 GHz) zenith is ~0.057 m.
  double x_band = IonosphereDelayThinShell(PI / 2, 10.0, 8.4e9, 1000.0);
  double expected = 40.3 * (10.0 * 1e16) / (8.4e9 * 8.4e9);  // vertical, mapping = 1 at zenith
  REQUIRE_THAT(x_band, WithinRel(expected, 1e-6));

  // The delay is dispersive (1/f^2): S-band (2.2 GHz) is (8.4/2.2)^2 larger than X-band.
  double s_band = IonosphereDelayThinShell(PI / 2, 10.0, 2.2e9, 1000.0);
  REQUIRE_THAT(s_band / x_band, WithinRel((8.4e9 * 8.4e9) / (2.2e9 * 2.2e9), 1e-6));

  // The slant delay grows toward the horizon (obliquity > 1).
  double low = IonosphereDelayThinShell(10.0 * RAD, 10.0, 8.4e9, 1000.0);
  REQUIRE(low > x_band);

  // Zero TEC or non-positive frequency yields no delay.
  REQUIRE(IonosphereDelayThinShell(PI / 2, 0.0, 8.4e9, 1000.0) == 0.0);
  REQUIRE(IonosphereDelayThinShell(PI / 2, 10.0, 0.0, 1000.0) == 0.0);
}

TEST_CASE("ground_station_corrections.shapiro") {
  // Earth-station (near geocenter) to lunar-distance receiver, along the Sun line.
  Vec3d earth(0.0, 0.0, 0.0);
  Vec3d sun(1.496e11, 0.0, 0.0);
  Vec3d tx(6.4e6, 0.0, 0.0);
  Vec3d rx(3.84e8, 0.0, 0.0);

  double earth_only = ShapiroRangeDelay(tx, rx, {{earth, GM_EARTH}});
  double with_sun = ShapiroRangeDelay(tx, rx, {{earth, GM_EARTH}, {sun, GM_SUN}});
  REQUIRE(earth_only > 0.0);
  REQUIRE(with_sun > earth_only);  // adding the Sun term increases the delay

  // Earth term matches the closed form 2 GM/c^2 * ln[(r1+r2+r12)/(r1+r2-r12)].
  double r1 = tx.norm(), r2 = rx.norm(), r12 = (rx - tx).norm();
  double expected = 2.0 * GM_EARTH / (C * C) * std::log((r1 + r2 + r12) / (r1 + r2 - r12));
  REQUIRE_THAT(earth_only, WithinRel(expected, 1e-9));

  // No bodies -> no delay.
  REQUIRE(ShapiroRangeDelay(tx, rx, {}) == 0.0);
}

TEST_CASE("ground_station_corrections.solid_earth_tide") {
  Vec3d station(R_EARTH, 0.0, 0.0);
  Vec3d moon(3.84e8, 0.0, 0.0);  // aligned with the station -> maximal radial bulge
  Vec3d sun(1.496e11, 0.0, 0.0);

  Vec3d disp
      = SolidEarthTideDisplacement(station, {{moon, GM_MOON}, {sun, GM_SUN}}, GM_EARTH, R_EARTH);
  // Sub-decimetre-to-decimetre scale, well under half a metre.
  REQUIRE(disp.norm() < 0.5);
  // Aligned geometry raises the crust radially outward (+x).
  REQUIRE(disp.x() > 0.05);

  // With the tide bodies at 90 deg the radial term flips sign (crust drawn inward).
  Vec3d moon_perp(0.0, 3.84e8, 0.0);
  Vec3d disp_perp = SolidEarthTideDisplacement(station, {{moon_perp, GM_MOON}}, GM_EARTH, R_EARTH);
  REQUIRE(disp_perp.x() < 0.0);

  // No tide-raising bodies -> no displacement.
  REQUIRE(SolidEarthTideDisplacement(station, {}, GM_EARTH, R_EARTH).norm() == 0.0);
}
