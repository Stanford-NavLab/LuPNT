#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"
using namespace lupnt;
using namespace Catch::Matchers;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

const Real M_ref = 0.06981317008;
const Real e_ref = 0.72000000000;
const Real E_ref = 0.24318719637;
const Real nu_ref = MeanToTrueAnomaly(M_ref, e_ref);

TEST_CASE("physics.orbit_state.anomaly") {
  const double eps = 1e-9;

  // Example 2-2 (Kepler’s equation)
  const Real M_ref = 0.06981317008;
  const Real e_ref = 0.72000000000;
  const Real E_ref = 0.24318719637;
  const Real nu_ref = MeanToTrueAnomaly(M_ref, e_ref);

  RequireNear(MeanToEccAnomaly(M_ref, e_ref), E_ref, eps);
  RequireNear(EccToMeanAnomaly(E_ref, e_ref), M_ref, eps);
  RequireNear(EccToTrueAnomaly(E_ref, e_ref), nu_ref, eps);
  RequireNear(TrueToEccAnomaly(nu_ref, e_ref), E_ref, eps);
  RequireNear(MeanToTrueAnomaly(M_ref, e_ref), nu_ref, eps);
  RequireNear(TrueToMeanAnomaly(nu_ref, e_ref), M_ref, eps);
}

TEST_CASE("physics.orbit_state.conversions") {
  const double eps = 1e-5;

  // Example 2-3 (Osculating Elements)
  Vec6 rv_ref(10e3, 40e3, -5e3, -1.5, 1, -0.1);  // [km, km/s]
  Vec6 coe_ref(25015.181022316396, 0.707977170662, 6.970729208731, 173.290163192243,
               91.552887356747, 144.224991174457);
  coe_ref.segment(2, 4) *= RAD;
  // [km, -, rad, rad, rad, rad]

  Real GM_earth_km = GM_EARTH * 1e-9;
  Vec6 rv = ClassicalToCart(coe_ref, GM_earth_km);
  RequireNear(rv, rv_ref, eps);

  Vec6 coe = CartToClassical(rv, GM_earth_km);
  RequireNear(coe(1), coe_ref(1), eps);
  RequireNear(coe, coe_ref, eps);
}

TEST_CASE("physics.orbit_state.conversions2") {
  double eps = 1e-2;

  Vec3 lla(41.386464272846, 2.169701993823, 123.456321000000003);
  Vec3 enu(3544.169276615858053, 8202.463864965198809, 9149.715157235310471);
  Vec3 aer = EastNorthUpToAzElRange(enu);
  Vec3 r_gs = LatLonAltToCart(lla, R_EARTH, WGS84_F);
  Vec3 r_sat = EastNorthUpToCart(enu, r_gs, R_EARTH, WGS84_F);

  RequireNear(aer, EastNorthUpToAzElRange(enu), eps);
  RequireNear(enu, AzElRangeToEastNorthUp(aer), eps);

  RequireNear(r_gs, LatLonAltToCart(lla, R_EARTH, WGS84_F), eps);
  Vec3 lla_roundtrip = CartToLatLonAlt(r_gs, R_EARTH, WGS84_F);
  RequireNear(lla(0), lla_roundtrip(0), eps);
  RequireNear(lla(1), lla_roundtrip(1), eps);

  RequireNear(r_sat, EastNorthUpToCart(enu, r_gs, R_EARTH, WGS84_F), eps);
  RequireNear(enu, CartToEastNorthUp(r_sat, r_gs, R_EARTH, WGS84_F), eps);

  RequireNear(aer, CartToAzElRange(r_sat, r_gs, R_EARTH, WGS84_F), eps);
  RequireNear(r_sat, AzElRangeToCart(aer, r_gs, R_EARTH, WGS84_F), eps);
}

TEST_CASE("physics.orbit_state.utils") {
  double eps = 1e-9;
  // Example 2-2 (Kepler’s equation)
  Real M, E;
  Real e = 0.72;

  M = 4 * RAD;
  E = MeanToEccAnomaly(M, e);
  RequireNear(E, 0.24318719629, 1e-10);

  M = 50 * RAD;
  E = MeanToEccAnomaly(M, e);
  RequireNear(E, 1.59249513093, 1e-10);

  // Example 2-3 (Osculating Elements)
  Vec6 rv(10e3, 40e3, -5e3, -1.5, 1, -0.1);  // [km, km/s]
  Real GM_earth_km = GM_EARTH * 1e-9;
  Vec6 coe = CartToClassical(rv, GM_earth_km);
  RequireNear(coe[0], 25015.181, 1e-3);
  RequireNear(coe[1], 0.7079772, 1e-7);
  RequireNear(coe[2] * DEG, 6.971, 1e-3);
  RequireNear(coe[3] * DEG, 173.290, 1e-3);
  RequireNear(coe[4] * DEG, 91.553, 1e-3);
  RequireNear(coe[5] * DEG, 144.225, 1e-3);

  // Example 2-4 (Topocentric satellite motion)
  EopFileData* eop_data = GetEopFileData();
  // 1962 1 1  37665  -0.012700   0.213000   0.0326338   0.0017230 0.064261
  // 0.006067 0.030000   0.030000  0.0020000  0.0014000    0.012000 0.002000
  RequireNear(eop_data->years(0), 1962, eps);
  RequireNear(eop_data->months(0), 1, eps);
  RequireNear(eop_data->days(0), 1, eps);
  RequireNear(eop_data->mjds_utc(0), 37665, eps);
  RequireNear(eop_data->x(0), -0.012700, eps);
  RequireNear(eop_data->y(0), 0.213000, eps);
  RequireNear(eop_data->ut1_utc(0), 0.0326338, eps);
  RequireNear(eop_data->lod(0), 0.0017230, eps);
  RequireNear(eop_data->dpsi(0), 0.064261, eps);
  RequireNear(eop_data->deps(0), 0.006067, eps);
  RequireNear(eop_data->xErr(0), 0.030000, eps);
  RequireNear(eop_data->yErr(0), 0.030000, eps);
  RequireNear(eop_data->ut1_utc_err(0), 0.0020000, eps);
  RequireNear(eop_data->lod_err(0), 0.0014000, eps);
  RequireNear(eop_data->dpsi_err(0), 0.012000, eps);
  RequireNear(eop_data->deps_err(0), 0.002000, eps);

  Real mjd_utc_1 = 37665;
  Real mjd_utc_2 = 37666;

  EopData res1;
  res1.x_pole = -0.012700 * RAD_ARCSEC;
  res1.y_pole = 0.213000 * RAD_ARCSEC;
  res1.ut1_utc = 0.0326338;
  res1.lod = 0.0017230;
  res1.dpsi = 0.064261 * RAD_ARCSEC;
  res1.deps = 0.006067 * RAD_ARCSEC;
  res1.sigma_x_pole = 0.030000 * RAD_ARCSEC;
  res1.sigma_y_pole = 0.030000 * RAD_ARCSEC;
  res1.sigma_ut1_utc = 0.0020000;

  EopData res2;
  res2.x_pole = -0.015900 * RAD_ARCSEC;
  res2.y_pole = 0.214100 * RAD_ARCSEC;
  res2.ut1_utc = 0.0320547;
  res2.lod = 0.0016690;
  res2.dpsi = 0.063979 * RAD_ARCSEC;
  res2.deps = 0.006290 * RAD_ARCSEC;
  res2.sigma_x_pole = 0.030000 * RAD_ARCSEC;
  res2.sigma_y_pole = 0.030000 * RAD_ARCSEC;
  res2.sigma_ut1_utc = 0.0020000;

  // Interpolation
  Real s = 0.4;
  Real mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
  eps = 1e-3;

  auto interp = [](Real x0, Real x1, Real s) { return x0 + (x1 - x0) * s; };
  EopData result = GetEopData(mjd_utc);
  RequireNear(result.x_pole, interp(res1.x_pole, res2.x_pole, s), eps);
  RequireNear(result.y_pole, interp(res1.y_pole, res2.y_pole, s), eps);
  RequireNear(result.ut1_utc, interp(res1.ut1_utc, res2.ut1_utc, s), eps);
  RequireNear(result.lod, interp(res1.lod, res2.lod, s), eps);
  RequireNear(result.dpsi, interp(res1.dpsi, res2.dpsi, s), eps);
  RequireNear(result.deps, interp(res1.deps, res2.deps, s), eps);
  RequireNear(result.sigma_x_pole, interp(res1.sigma_x_pole, res2.sigma_x_pole, s), eps);
  RequireNear(result.sigma_y_pole, interp(res1.sigma_y_pole, res2.sigma_y_pole, s), eps);
  RequireNear(result.sigma_ut1_utc, interp(res1.sigma_ut1_utc, res2.sigma_ut1_utc, s), eps);

  s = 0.2;
  mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
  result = GetEopData(mjd_utc);
  RequireNear(result.x_pole, res1.x_pole, eps);
  RequireNear(result.y_pole, res1.y_pole, eps);
  RequireNear(result.ut1_utc, res1.ut1_utc, eps);
  RequireNear(result.lod, res1.lod, eps);
  RequireNear(result.dpsi, res1.dpsi, eps);
  RequireNear(result.deps, res1.deps, eps);
  RequireNear(result.sigma_x_pole, res1.sigma_x_pole, eps);
  RequireNear(result.sigma_y_pole, res1.sigma_y_pole, eps);
  RequireNear(result.sigma_ut1_utc, res1.sigma_ut1_utc, eps);

  s = 0.8;
  mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
  result = GetEopData(mjd_utc);
  RequireNear(result.x_pole, res2.x_pole, eps);
  RequireNear(result.y_pole, res2.y_pole, eps);
  RequireNear(result.ut1_utc, res2.ut1_utc, eps);
  RequireNear(result.lod, res2.lod, eps);
  RequireNear(result.dpsi, res2.dpsi, eps);
  RequireNear(result.deps, res2.deps, eps);
  RequireNear(result.sigma_x_pole, res2.sigma_x_pole, eps);
  RequireNear(result.sigma_y_pole, res2.sigma_y_pole, eps);
  RequireNear(result.sigma_ut1_utc, res2.sigma_ut1_utc, eps);

  // Ground station
  Real lon_gs = 11;
  Real lat_gs = 48;
  Real alt_gs = 0;
  Vec3 r_gs = LatLonAltToCart(Vec3(lat_gs, lon_gs, alt_gs), R_EARTH, WGS84_F);
  Vec3 r_gs_ref(4197.16082495916e3, 815.845418656284e3, 4716.87633011541e3);
  RequireNear(r_gs, r_gs_ref, eps);
  Vec3 r_geod = CartToLatLonAlt(r_gs, R_EARTH, WGS84_F);
  RequireNear(r_geod(0), 48, 1e-2);
  RequireNear(r_geod(1), 11, 1e-2);

  // Spacecraft
  Real a = 960e3 + R_EARTH;              // Semimajor axis [m]
  e = 0;                                 // Eccentricity
  Real i = 97 * RAD;                     // Inclination [rad]
  Real Omega = 130.7 * RAD;              // RA ascend. node [rad]
  Real omega = 0 * RAD;                  // Argument of latitude [rad]
  Real M0 = 0 * RAD;                     // Mean anomaly at epoch [rad]
  Vec6 coe0(a, e, i, Omega, omega, M0);  // Classical orbital elements
  Vec6 coe0_ref(7338.137e3, 0, 1.6929693744345, 2.28114533235659, 0, 0);
  RequireNear(coe0, coe0_ref, eps);

  // // Propagation
  // Real minute = 6;
  // mjd_utc = mjd0_utc + minute * SECS_MINUTE / SECS_DAY;
  // Real mjd_ut1 = UtcToUt1(mjd_utc);
  // Real t = (mjd_utc - mjd0_utc) * SECS_DAY;
  // coe = KeplerianDynamics::PropagateClassicalOE(coe0, t, GM_EARTH);
  // Vec6 rv_eci = ClassicalToCart(coe, GM_EARTH);
  // Real theta_era = GreenwichMeanSideRealTime(mjd_ut1);
  // Mat3 eci_to_ecef = RotZ(theta_era);
  // Vec3 r_ecef = eci_to_ecef * rv_eci.head(3);
  // auto [az, el, rho] = unpack(CartToAzElRange(r_gs, r_ecef, WGS84_F));

  // Real mjd_ut1_ref = 50449.0041653811;
  // Vec6 rv_eci_ref(-4235.95304225382, 5409.87676038481, 2576.46849234333, 2.33703511947809,
  //                 -1.42872140442756, 6.84222523949854);
  // Vec3 r_ecef_ref(6182.21120035545, 2998.38780229154, 2576.46849234333);
  // Real az_ref = 2.63651291242072;
  // Real el_ref = -0.00222785822255024;
  // Real rho_ref = 3644.89532925451;
  // RequireNear(mjd_ut1, mjd_ut1_ref, eps);
  // RequireNear(rv_eci, rv_eci_ref, eps);
  // RequireNear(r_ecef, r_ecef_ref, eps);
  // RequireNear(az, az_ref, eps);
  // RequireNear(el, el_ref, eps);
  // RequireNear(rho, rho_ref, eps);
}
