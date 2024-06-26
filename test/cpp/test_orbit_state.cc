#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "test_utils.cc"
using namespace lupnt;
using namespace Catch::Matchers;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

TEST_CASE("OrbitState", "Utils") {
  double eps = 1e-9;
  // Example 2-2 (Kepler’s equation)
  real M, E;
  real e = 0.72;

  M = 4 * RAD_PER_DEG;
  E = MeanToEccentricAnomaly(M, e);
  RequireNearReal(E, 0.24318719629, 1e-10);

  M = 50 * RAD_PER_DEG;
  E = MeanToEccentricAnomaly(M, e);
  RequireNearReal(E, 1.59249513093, 1e-10);

  // Example 2-3 (Osculating Elements)
  Vector6 rv(10e3, 40e3, -5e3, -1.5, 1, -0.1);  // [km, km/s]
  Vector6 coe = CartesianToClassical(rv, GM_EARTH);
  RequireNearReal(coe[0], 25015.181, 1e-3);
  RequireNearReal(coe[1], 0.7079772, 1e-7);
  RequireNearReal(coe[2] * DEG_PER_RAD, 6.971, 1e-3);
  RequireNearReal(coe[3] * DEG_PER_RAD, 173.290, 1e-3);
  RequireNearReal(coe[4] * DEG_PER_RAD, 91.553, 1e-3);
  RequireNearReal(coe[5] * DEG_PER_RAD, 144.225, 1e-3);

  // Example 2-4 (Topocentric satellite motion)
  std::shared_ptr<EOPData> eop_data =
      LoadEOPData(GetFilePath("eopc04_08.62-now"));
  // 1962 1 1  37665  -0.012700   0.213000   0.0326338   0.0017230 0.064261
  // 0.006067 0.030000   0.030000  0.0020000  0.0014000    0.012000 0.002000
  RequireNearReal(eop_data->years(0), 1962, eps);
  RequireNearReal(eop_data->months(0), 1, eps);
  RequireNearReal(eop_data->days(0), 1, eps);
  RequireNearReal(eop_data->mjds(0), 37665, eps);
  RequireNearReal(eop_data->x(0), -0.012700, eps);
  RequireNearReal(eop_data->y(0), 0.213000, eps);
  RequireNearReal(eop_data->ut1_utc(0), 0.0326338, eps);
  RequireNearReal(eop_data->lod(0), 0.0017230, eps);
  RequireNearReal(eop_data->dPsi(0), 0.064261, eps);
  RequireNearReal(eop_data->dEps(0), 0.006067, eps);
  RequireNearReal(eop_data->xErr(0), 0.030000, eps);
  RequireNearReal(eop_data->yErr(0), 0.030000, eps);
  RequireNearReal(eop_data->ut1_utcErr(0), 0.0020000, eps);
  RequireNearReal(eop_data->lodErr(0), 0.0014000, eps);
  RequireNearReal(eop_data->dPsiErr(0), 0.012000, eps);
  RequireNearReal(eop_data->dEpsErr(0), 0.002000, eps);

  real mjd_utc_1 = 37665;
  real mjd_utc_2 = 37666;

  EOPResult res1;
  res1.x_pole = -0.012700 * RAD_PER_ARCSEC;
  res1.y_pole = 0.213000 * RAD_PER_ARCSEC;
  res1.UT1_UTC = 0.0326338;
  res1.LOD = 0.0017230;
  res1.dPsi = 0.064261 * RAD_PER_ARCSEC;
  res1.dEps = 0.006067 * RAD_PER_ARCSEC;
  res1.dx_pole = 0.030000 * RAD_PER_ARCSEC;
  res1.dy_pole = 0.030000 * RAD_PER_ARCSEC;
  res1.TAI_UTC = 0.0020000;

  EOPResult res2;
  res2.x_pole = -0.015900 * RAD_PER_ARCSEC;
  res2.y_pole = 0.214100 * RAD_PER_ARCSEC;
  res2.UT1_UTC = 0.0320547;
  res2.LOD = 0.0016690;
  res2.dPsi = 0.063979 * RAD_PER_ARCSEC;
  res2.dEps = 0.006290 * RAD_PER_ARCSEC;
  res2.dx_pole = 0.030000 * RAD_PER_ARCSEC;
  res2.dy_pole = 0.030000 * RAD_PER_ARCSEC;
  res2.TAI_UTC = 0.0020000;

  // Interpolation
  real s = 0.4;
  real mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
  eps = 1e-6;

  auto interp = [](real x0, real x1, real s) { return x0 + (x1 - x0) * s; };
  EOPResult result = InterpolateEOPData(eop_data, mjd_utc, true);
  RequireNearReal(result.x_pole, interp(res1.x_pole, res2.x_pole, s), eps);
  RequireNearReal(result.y_pole, interp(res1.y_pole, res2.y_pole, s), eps);
  RequireNearReal(result.UT1_UTC, interp(res1.UT1_UTC, res2.UT1_UTC, s), eps);
  RequireNearReal(result.LOD, interp(res1.LOD, res2.LOD, s), eps);
  RequireNearReal(result.dPsi, interp(res1.dPsi, res2.dPsi, s), eps);
  RequireNearReal(result.dEps, interp(res1.dEps, res2.dEps, s), eps);
  RequireNearReal(result.dx_pole, interp(res1.dx_pole, res2.dx_pole, s), eps);
  RequireNearReal(result.dy_pole, interp(res1.dy_pole, res2.dy_pole, s), eps);
  RequireNearReal(result.TAI_UTC, interp(res1.TAI_UTC, res2.TAI_UTC, s), eps);

  s = 0.2;
  mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
  result = InterpolateEOPData(eop_data, mjd_utc, false);
  RequireNearReal(result.x_pole, res1.x_pole, eps);
  RequireNearReal(result.y_pole, res1.y_pole, eps);
  RequireNearReal(result.UT1_UTC, res1.UT1_UTC, eps);
  RequireNearReal(result.LOD, res1.LOD, eps);
  RequireNearReal(result.dPsi, res1.dPsi, eps);
  RequireNearReal(result.dEps, res1.dEps, eps);
  RequireNearReal(result.dx_pole, res1.dx_pole, eps);
  RequireNearReal(result.dy_pole, res1.dy_pole, eps);
  RequireNearReal(result.TAI_UTC, res1.TAI_UTC, eps);

  s = 0.8;
  mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
  result = InterpolateEOPData(eop_data, mjd_utc, false);
  RequireNearReal(result.x_pole, res2.x_pole, eps);
  RequireNearReal(result.y_pole, res2.y_pole, eps);
  RequireNearReal(result.UT1_UTC, res2.UT1_UTC, eps);
  RequireNearReal(result.LOD, res2.LOD, eps);
  RequireNearReal(result.dPsi, res2.dPsi, eps);
  RequireNearReal(result.dEps, res2.dEps, eps);
  RequireNearReal(result.dx_pole, res2.dx_pole, eps);
  RequireNearReal(result.dy_pole, res2.dy_pole, eps);
  RequireNearReal(result.TAI_UTC, res2.TAI_UTC, eps);

  // Reloading the EOP data points to the same instance
  std::shared_ptr<EOPData> eop_data_reloaded =
      LoadEOPData(GetFilePath("eopc04_08.62-now"));
  REQUIRE(eop_data == eop_data_reloaded);

  // Time
  real mjd0_utc = DateToModifiedJulianDate(1997, 1, 1, 0, 0, 0);
  RequireNearReal(mjd0_utc, 50449, eps);

  // Ground station
  real lon_gs = 11 * RAD_PER_DEG;
  real lat_gs = 48 * RAD_PER_DEG;
  real alt_gs = 0;
  Vector3 r_gs = GeodeticToCartesian(Vector3(lat_gs, lon_gs, alt_gs), R_EARTH,
                                     FLATTENING_EARTH_WGS84);
  Vector3 r_gs_ref(4197.16082495916, 815.845418656284, 4716.87633011541);
  RequireNearRealVec(r_gs, r_gs_ref, eps);
  Vector3 r_geod = CartesianToGeodetic(r_gs, R_EARTH, FLATTENING_EARTH_WGS84);
  Vector3 r_geod_ref(48 * RAD_PER_DEG, 11 * RAD_PER_DEG, 0);
  RequireNearRealVec(r_geod, r_geod_ref, eps);

  // Spacecraft
  real a = 960 + R_EARTH;                   // Semimajor axis [Km]
  e = 0;                                    // Eccentricity
  real i = 97 * RAD_PER_DEG;                // Inclination [rad]
  real Omega = 130.7 * RAD_PER_DEG;         // RA ascend. node [rad]
  real omega = 0 * RAD_PER_DEG;             // Argument of latitude [rad]
  real M0 = 0 * RAD_PER_DEG;                // Mean anomaly at epoch [rad]
  Vector6 coe0(a, e, i, Omega, omega, M0);  // Classical orbital elements
  Vector6 coe0_ref(7338.137, 0, 1.6929693744345, 2.28114533235659, 0, 0);
  RequireNearRealVec(coe0, coe0_ref, eps);

  // Propagation
  real minute = 6;
  mjd_utc = mjd0_utc + minute * SECS_PER_MINUTE / SECS_PER_DAY;
  real mjd_ut1 = UTCtoUT1(mjd_utc);
  real t = (mjd_utc - mjd0_utc) * SECS_PER_DAY;
  coe = KeplerianDynamics::PropagateClassicalOE(coe0, t, GM_EARTH);
  Vector6 rv_eci = ClassicalToCartesian(coe, GM_EARTH);
  real theta_era = GreenwichMeanSiderealTime(mjd_ut1);
  Matrix3 eci_to_ecef = Rot3(theta_era);
  Vector3 r_ecef = eci_to_ecef * rv_eci.head(3);
  auto [az, el, rho] = unpack(
      CartesianToAzimuthElevationRange(r_gs, r_ecef, FLATTENING_EARTH_WGS84));

  real mjd_ut1_ref = 50449.0041653811;
  Vector6 rv_eci_ref(-4235.95304225382, 5409.87676038481, 2576.46849234333,
                     2.33703511947809, -1.42872140442756, 6.84222523949854);
  Vector3 r_ecef_ref(6182.21120035545, 2998.38780229154, 2576.46849234333);
  real az_ref = 2.63651291242072;
  real el_ref = -0.00222785822255024;
  real rho_ref = 3644.89532925451;
  RequireNearReal(mjd_ut1, mjd_ut1_ref, eps);
  RequireNearRealVec(rv_eci, rv_eci_ref, eps);
  RequireNearRealVec(r_ecef, r_ecef_ref, eps);
  RequireNearReal(az, az_ref, eps);
  RequireNearReal(el, el_ref, eps);
  RequireNearReal(rho, rho_ref, eps);
}