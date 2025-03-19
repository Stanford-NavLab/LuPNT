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
const Real nu_ref = Mean2TrueAnomaly(M_ref, e_ref);

TEST_CASE("Anomaly", "[OrbitState]") {
  const double eps = 1e-9;

  // Example 2-2 (Kepler’s equation)
  const Real M_ref = 0.06981317008;
  const Real e_ref = 0.72000000000;
  const Real E_ref = 0.24318719637;
  const Real nu_ref = Mean2TrueAnomaly(M_ref, e_ref);

  RequireNear(Mean2EccAnomaly(M_ref, e_ref), E_ref, eps);
  RequireNear(Ecc2MeanAnomaly(E_ref, e_ref), M_ref, eps);
  RequireNear(Ecc2TrueAnomaly(E_ref, e_ref), nu_ref, eps);
  RequireNear(True2EccAnomaly(nu_ref, e_ref), E_ref, eps);
  RequireNear(Mean2TrueAnomaly(M_ref, e_ref), nu_ref, eps);
  RequireNear(True2MeanAnomaly(nu_ref, e_ref), M_ref, eps);
}

TEST_CASE("Conversions", "[OrbitState]") {
  const double eps = 1e-5;

  // Example 2-3 (Osculating Elements)
  Vec6 rv_ref(10e3, 40e3, -5e3, -1.5, 1, -0.1);  // [km, km/s]
  Vec6 coe_ref(25015.181022316396, 0.707977170662, 6.970729208731, 173.290163192243,
               91.552887356747, 144.224991174457);
  coe_ref.segment(2, 4) *= RAD;
  // [km, -, rad, rad, rad, rad]

  Vec6 rv = Classical2Cart(coe_ref, GM_EARTH);
  RequireNear(rv, rv_ref, eps);

  Vec6 coe = Cart2Classical(rv, GM_EARTH);
  RequireNear(coe(1), coe_ref(1), eps);
  RequireNear(coe, coe_ref, eps);
}

TEST_CASE("OrbitState", "Coordinates") {
  double eps = 1e-8;

  // Vec3 lla(0.722346906730713, 0.037849373172501, 123.456321000000003);
  // Vec3 aer(0.407857156994684, 0.797247419062867, 12789.012345600000117);
  // Vec3 enu(3544.169276615858053, 8202.463864965198809, 9149.715157235310471);
  // Vec3 r_gs(4881.378514959144923, 184.845393739831565, 4276.423188333580583);
  // Vec3 r_sat(6187.826901711936443, 3781.026709647304870, 16479.681332443218707);

  // RequireNear(aer, EastNorthUp2AzElRange(enu), eps);
  // RequireNear(enu, AzElRange2EastNorthUp(aer), eps);

  // RequireNear(r_gs, LatLonAlt2Cart(lla, R_EARTH, WGS84_F), eps);
  // RequireNear(lla, Cart2LatLonAlt(r_gs, R_EARTH, WGS84_F), eps);

  // RequireNear(r_sat, EastNorthUp2Cart(enu, r_gs, R_EARTH, WGS84_F), eps);
  // RequireNear(enu, Cart2EastNorthUp(r_sat, r_gs, R_EARTH, WGS84_F), eps);

  // RequireNear(aer, Cart2AzElRange(r_sat, r_gs, R_EARTH, WGS84_F), eps);
  // RequireNear(r_sat, AzElRange2Cart(aer, r_gs, R_EARTH, WGS84_F), eps);
}

// TEST_CASE("OrbitState", "Utils") {
//   double eps = 1e-9;
//   // Example 2-2 (Kepler’s equation)
//   real M, E;
//   real e = 0.72;

//   M = 4 * RAD_PER_DEG;
//   E = Mean2EccAnomaly(M, e);
//   RequireNear(E, 0.24318719629, 1e-10);

//   M = 50 * RAD_PER_DEG;
//   E = Mean2EccAnomaly(M, e);
//   RequireNear(E, 1.59249513093, 1e-10);

//   // Example 2-3 (Osculating Elements)
//   Vec6 rv(10e3, 40e3, -5e3, -1.5, 1, -0.1);  // [km, km/s]
//   Vec6 coe = Cart2Classical(rv, GM_EARTH);
//   RequireNear(coe[0], 25015.181, 1e-3);
//   RequireNear(coe[1], 0.7079772, 1e-7);
//   RequireNear(coe[2] * DEG_PER_RAD, 6.971, 1e-3);
//   RequireNear(coe[3] * DEG_PER_RAD, 173.290, 1e-3);
//   RequireNear(coe[4] * DEG_PER_RAD, 91.553, 1e-3);
//   RequireNear(coe[5] * DEG_PER_RAD, 144.225, 1e-3);

//   // Example 2-4 (Topocentric satellite motion)
//   Ptr<EopFileData> eop_data =
//       LoadEopFileData(GetFilePath("eopc04_08.62-now"));
//   // 1962 1 1  37665  -0.012700   0.213000   0.0326338   0.0017230 0.064261
//   // 0.006067 0.030000   0.030000  0.0020000  0.0014000    0.012000 0.002000
//   RequireNear(eop_data->years(0), 1962, eps);
//   RequireNear(eop_data->months(0), 1, eps);
//   RequireNear(eop_data->days(0), 1, eps);
//   RequireNear(eop_data->mjds_utc(0), 37665, eps);
//   RequireNear(eop_data->x(0), -0.012700, eps);
//   RequireNear(eop_data->y(0), 0.213000, eps);
//   RequireNear(eop_data->ut1_utc(0), 0.0326338, eps);
//   RequireNear(eop_data->lod(0), 0.0017230, eps);
//   RequireNear(eop_data->dpsi(0), 0.064261, eps);
//   RequireNear(eop_data->deps(0), 0.006067, eps);
//   RequireNear(eop_data->xErr(0), 0.030000, eps);
//   RequireNear(eop_data->yErr(0), 0.030000, eps);
//   RequireNear(eop_data->ut1_utc_err(0), 0.0020000, eps);
//   RequireNear(eop_data->lod_err(0), 0.0014000, eps);
//   RequireNear(eop_data->dpsi_err(0), 0.012000, eps);
//   RequireNear(eop_data->deps_err(0), 0.002000, eps);

//   real mjd_utc_1 = 37665;
//   real mjd_utc_2 = 37666;

//   EopData res1;
//   res1.x_pole = -0.012700 * RAD_PER_ARCSEC;
//   res1.y_pole = 0.213000 * RAD_PER_ARCSEC;
//   res1.ut1_utc = 0.0326338;
//   res1.lod = 0.0017230;
//   res1.dpsi = 0.064261 * RAD_PER_ARCSEC;
//   res1.deps = 0.006067 * RAD_PER_ARCSEC;
//   res1.dx_pole = 0.030000 * RAD_PER_ARCSEC;
//   res1.dy_pole = 0.030000 * RAD_PER_ARCSEC;
//   res1.tai_utc = 0.0020000;

//   EopData res2;
//   res2.x_pole = -0.015900 * RAD_PER_ARCSEC;
//   res2.y_pole = 0.214100 * RAD_PER_ARCSEC;
//   res2.ut1_utc = 0.0320547;
//   res2.lod = 0.0016690;
//   res2.dpsi = 0.063979 * RAD_PER_ARCSEC;
//   res2.deps = 0.006290 * RAD_PER_ARCSEC;
//   res2.dx_pole = 0.030000 * RAD_PER_ARCSEC;
//   res2.dy_pole = 0.030000 * RAD_PER_ARCSEC;
//   res2.tai_utc = 0.0020000;

//   // Interpolation
//   real s = 0.4;
//   real mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
//   eps = 1e-6;

//   auto interp = [](real x0, real x1, real s) { return x0 + (x1 - x0) * s; };
//   EopData result = GetEopData(eop_data, mjd_utc, true);
//   RequireNear(result.x_pole, interp(res1.x_pole, res2.x_pole, s), eps);
//   RequireNear(result.y_pole, interp(res1.y_pole, res2.y_pole, s), eps);
//   RequireNear(result.ut1_utc, interp(res1.ut1_utc, res2.ut1_utc, s),
//   eps); RequireNear(result.lod, interp(res1.lod, res2.lod, s), eps);
//   RequireNear(result.dpsi, interp(res1.dpsi, res2.dpsi, s), eps);
//   RequireNear(result.deps, interp(res1.deps, res2.deps, s), eps);
//   RequireNear(result.dx_pole, interp(res1.dx_pole, res2.dx_pole, s),
//   eps); RequireNear(result.dy_pole, interp(res1.dy_pole, res2.dy_pole,
//   s), eps); RequireNear(result.tai_utc, interp(res1.tai_utc,
//   res2.tai_utc, s), eps);

//   s = 0.2;
//   mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
//   result = GetEopData(eop_data, mjd_utc, false);
//   RequireNear(result.x_pole, res1.x_pole, eps);
//   RequireNear(result.y_pole, res1.y_pole, eps);
//   RequireNear(result.ut1_utc, res1.ut1_utc, eps);
//   RequireNear(result.lod, res1.lod, eps);
//   RequireNear(result.dpsi, res1.dpsi, eps);
//   RequireNear(result.deps, res1.deps, eps);
//   RequireNear(result.dx_pole, res1.dx_pole, eps);
//   RequireNear(result.dy_pole, res1.dy_pole, eps);
//   RequireNear(result.tai_utc, res1.tai_utc, eps);

//   s = 0.8;
//   mjd_utc = mjd_utc_1 + s * (mjd_utc_2 - mjd_utc_1);
//   result = GetEopData(eop_data, mjd_utc, false);
//   RequireNear(result.x_pole, res2.x_pole, eps);
//   RequireNear(result.y_pole, res2.y_pole, eps);
//   RequireNear(result.ut1_utc, res2.ut1_utc, eps);
//   RequireNear(result.lod, res2.lod, eps);
//   RequireNear(result.dpsi, res2.dpsi, eps);
//   RequireNear(result.deps, res2.deps, eps);
//   RequireNear(result.dx_pole, res2.dx_pole, eps);
//   RequireNear(result.dy_pole, res2.dy_pole, eps);
//   RequireNear(result.tai_utc, res2.tai_utc, eps);

//   // Reloading the EOP data points to the same instance
//   Ptr<EopFileData> eop_data_reloaded =
//       LoadEopFileData(GetFilePath("eopc04_08.62-now"));
//   REQUIRE(eop_data == eop_data_reloaded);

//   // Time
//   real mjd0_utc = Date2ModifiedJulianDate(1997, 1, 1, 0, 0, 0);
//   RequireNear(mjd0_utc, 50449, eps);

//   // Ground station
//   real lon_gs = 11 * RAD_PER_DEG;
//   real lat_gs = 48 * RAD_PER_DEG;
//   real alt_gs = 0;
//   Vec3 r_gs = Geodetic2Cart(Vec3(lat_gs, lon_gs, alt_gs), R_EARTH,
//                             FLATTENING_EARTH_WGS84);
//   Vec3 r_gs_ref(4197.16082495916, 815.845418656284, 4716.87633011541);
//   RequireNear(r_gs, r_gs_ref, eps);
//   Vec3 r_geod = Cart2Geodetic(r_gs, R_EARTH, FLATTENING_EARTH_WGS84);
//   Vec3 r_geod_ref(48 * RAD_PER_DEG, 11 * RAD_PER_DEG, 0);
//   RequireNear(r_geod, r_geod_ref, eps);

//   // Spacecraft
//   real a = 960 + R_EARTH;                // Semimajor axis [Km]
//   e = 0;                                 // Eccentricity
//   real i = 97 * RAD_PER_DEG;             // Inclination [rad]
//   real Omega = 130.7 * RAD_PER_DEG;      // RA ascend. node [rad]
//   real omega = 0 * RAD_PER_DEG;          // Argument of latitude [rad]
//   real M0 = 0 * RAD_PER_DEG;             // Mean anomaly at epoch [rad]
//   Vec6 coe0(a, e, i, Omega, omega, M0);  // Classical orbital elements
//   Vec6 coe0_ref(7338.137, 0, 1.6929693744345, 2.28114533235659, 0, 0);
//   RequireNear(coe0, coe0_ref, eps);

//   // Propagation
//   real minute = 6;
//   mjd_utc = mjd0_utc + minute * SECS_MINUTE / SECS_PER_DAY;
//   real mjd_ut1 = UTC2UT1(mjd_utc);
//   real t = (mjd_utc - mjd0_utc) * SECS_PER_DAY;
//   coe = KeplerianDynamics::PropagateClassicalOE(coe0, t, GM_EARTH);
//   Vec6 rv_eci = Classical2Cart(coe, GM_EARTH);
//   real theta_era = GreenwichMeanSiderealTime(mjd_ut1);
//   Mat3 eci_to_ecef = RotZ(theta_era);
//   Vec3 r_ecef = eci_to_ecef * rv_eci.head(3);
//   auto [az, el, rho] =
//       unpack(Cart2AzElRange(r_gs, r_ecef, FLATTENING_EARTH_WGS84));

//   real mjd_ut1_ref = 50449.0041653811;
//   Vec6 rv_eci_ref(-4235.95304225382, 5409.87676038481, 2576.46849234333,
//                   2.33703511947809, -1.42872140442756, 6.84222523949854);
//   Vec3 r_ecef_ref(6182.21120035545, 2998.38780229154, 2576.46849234333);
//   real az_ref = 2.63651291242072;
//   real el_ref = -0.00222785822255024;
//   real rho_ref = 3644.89532925451;
//   RequireNear(mjd_ut1, mjd_ut1_ref, eps);
//   RequireNear(rv_eci, rv_eci_ref, eps);
//   RequireNear(r_ecef, r_ecef_ref, eps);
//   RequireNear(az, az_ref, eps);
//   RequireNear(el, el_ref, eps);
//   RequireNear(rho, rho_ref, eps);
// }
