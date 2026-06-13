#include "lupnt/conversions/coordinate_conversions.h"

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/states/state.h"
namespace lupnt {

  /**
   * @brief Convert spherical coordinates to Cartesian coordinates
   * @param enu East-North-Up coordinates (e, n, u) [m]
   * @return Azimuth-Elevation-Range coordinates (az, el, rho) [rad, rad, m]
   */
  State EastNorthUpToAzElRange(const State& enu) {
    Real e = enu(0), n = enu(1), u = enu(2);
    Real range = enu.norm();
    Real az = atan2(e, n);
    Real el = asin(u / range);
    return AzElRange({az, el, range}, enu.GetFrame());
  }
  VEC_IMP_VECTOR(EastNorthUpToAzElRange, 3);

  /**
   * @brief Convert Cartesian coordinates to spherical coordinates
   * @param aer Azimuth-Elevation-Range coordinates (az, el, rho) [rad, rad, m]
   * @return East-North-Up coordinates (e, n, u) [m]
   */
  State AzElRangeToEastNorthUp(const State& aer) {
    Real az = aer(0), el = aer(1), range = aer(2);
    Real e = range * cos(el) * sin(az);
    Real n = range * cos(el) * cos(az);
    Real u = range * sin(el);
    return Cart3({e, n, u}, aer.GetFrame());
  }
  VEC_IMP_VECTOR(AzElRangeToEastNorthUp, 3);

  /**
   * @brief Convert latitude, longitude, and altitude to Cartesian coordinates
   * @param lla Latitude, longitude, and altitude (lat, lon, alt) [deg, deg, m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Cartesian coordinates (x, y, z) [m]
   */
  State LatLonAltToCart(const State& lla, Real R_body, Real flattening) {
    Real lat = lla(0) * RAD, lon = lla(1) * RAD, alt = lla(2);
    Real e2 = flattening * (2 - flattening);
    Real cos_lat = cos(lat);
    Real sin_lat = sin(lat);
    Real N = R_body / sqrt(1 - e2 * sin_lat * sin_lat);
    Real x = (N + alt) * cos_lat * cos(lon);
    Real y = (N + alt) * cos_lat * sin(lon);
    Real z = ((1 - e2) * N + alt) * sin_lat;
    return Cart3({x, y, z}, lla.GetFrame());
  }
  VEC_IMP_VECTOR_REAL(LatLonAltToCart, 3);

  /**
   * @brief Convert Cartesian coordinates to latitude, longitude, and altitude
   * @param cart Cartesian coordinates (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Latitude, longitude, and altitude [deg, deg, m]
   */
  State CartToLatLonAlt(const State& cart, Real R_body, Real flattening) {
    if (abs(R_body) < EPS) R_body = cart.norm();

    Real x = cart(0), y = cart(1), z = cart(2);
    Real e2 = flattening * (2 - flattening);
    const int max_iterations = 1000;

    Real delta_z = e2 * z;  // Initial guess for Δz
    Real delta_z_ant = delta_z + 1e3 * EPS;
    Real sin_phi, N;
    int it = 0;
    while (it < max_iterations && abs(delta_z_ant - delta_z) >= EPS) {
      delta_z_ant = delta_z;
      sin_phi = (z + delta_z) / sqrt(x * x + y * y + (z + delta_z) * (z + delta_z));
      N = R_body / sqrt(1 - e2 * sin_phi * sin_phi);
      delta_z = N * e2 * sin_phi;
      it++;
    }
    if (it >= max_iterations) {
      LUPNT_CHECK(false, "CartToLatLonAlt did not converge", "CartToLatLonAlt");
    }

    Real lon = atan2(y, x);
    Real lat = atan2(z + delta_z, sqrt(x * x + y * y));
    Real h = sqrt(x * x + y * y + (z + delta_z) * (z + delta_z)) - N;
    return LatLonAlt({lat * DEG, lon * DEG, h}, cart.GetFrame());
  }
  VEC_IMP_VECTOR_REAL(CartToLatLonAlt, 3);

  /**
   * @brief Convert East-North-Up to Cartesian rotation matrix
   * @param xyz_ref Reference position vector (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Rotation matrix from East-North-Up to Cartesian
   */
  Mat3 RotEastNorthUpToCart(const State& xyz_ref, Real R_body, Real flattening) {
    return RotCartToEastNorthUp(xyz_ref, R_body, flattening).transpose();
  }

  /**
   * @brief Convert Cartesian to East-North-Up rotation matrix
   * @param xyz_ref Reference position vector (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Rotation matrix from Cartesian to East-North-Up
   */
  Mat3 RotCartToEastNorthUp(const State& xyz_ref, Real R_body, Real flattening) {
    Vec3 lla = CartToLatLonAlt(xyz_ref, R_body, flattening);
    Real lat = lla(0) * RAD, lon = lla(1) * RAD;
    Mat3 R_body2enu = RotX(PI_OVER_TWO - lat) * RotZ(PI_OVER_TWO + lon);
    return R_body2enu;
  }

  /**
   * @brief Convert East-North-Up coordinates to Cartesian coordinates
   * @param enu East-North-Up coordinates (e, n, u) [m] or (e, n, u, ve, vn, vu) [m, m/s]
   * @param xyz_ref Reference position vector (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Cartesian coordinates (x, y, z) [m] or (x, y, z, vx, vy, vz) [m, m/s]
   */
  State EastNorthUpToCart(const State& enu, const State& xyz_ref, Real R_body, Real flattening) {
    Mat3 R_enu2body = RotEastNorthUpToCart(xyz_ref, R_body, flattening);
    Vec3 r = R_enu2body * enu + xyz_ref;
    if (enu.size() == 6) {
      Vec3 v = R_enu2body * enu.tail(3);
      return Cart6(r, v, enu.GetFrame());
    } else {
      return Cart3(r, enu.GetFrame());
    }
  }
  VEC_IMP_VECTOR_VECTOR(EastNorthUpToCart, 3);

  /**
   * @brief Convert Cartesian coordinates to East-North-Up coordinates
   * @param xyz Cartesian position vector (x, y, z) [m] or (x, y, z, vx, vy, vz) [m, m/s]
   * @param xyz_ref Reference position vector (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return East-North-Up coordinates (e, n, u) [m] or (e, n, u, ve, vn, vu) [m, m/s]
   */
  State CartToEastNorthUp(const State& xyz, const State& xyz_ref, Real R_body, Real flattening) {
    Mat3 R_body2enu = RotCartToEastNorthUp(xyz_ref, R_body, flattening);
    Vec3 r_enu = R_body2enu * (xyz - xyz_ref);
    if (xyz.size() == 6) {
      Vec3 v_enu = R_body2enu * xyz.tail(3);
      return Cart6(r_enu, v_enu, xyz.GetFrame());
    } else {
      return Cart3(r_enu, xyz.GetFrame());
    }
  }
  VEC_IMP_VECTOR_VECTOR(CartToEastNorthUp, 3);

  /**
   * @brief Convert Cartesian coordinates to Azimuth-Elevation-Range coordinates
   * @param xyz Cartesian position vector (x, y, z) [m]
   * @param xyz_ref Reference position vector (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Azimuth-Elevation-Range coordinates (az, el, rho) [rad, rad, m]
   */
  State CartToAzElRange(const State& xyz, const State& xyz_ref, Real R_body, Real flattening) {
    State enu = CartToEastNorthUp(xyz, xyz_ref, R_body, flattening);
    State aer = EastNorthUpToAzElRange(enu);
    return AzElRange(aer, xyz.GetFrame());
  }
  VEC_IMP_VECTOR_VECTOR(CartToAzElRange, 3);

  /**
   * @brief Convert Azimuth-Elevation-Range coordinates to Cartesian coordinates
   * @param aer Azimuth-Elevation-Range coordinates (az, el, rho) [rad, rad, m]
   * @param xyz_ref Reference position vector (x, y, z) [m]
   * @param R_body Body radius [m]
   * @param flattening Body flattening [-]
   * @return Cartesian coordinates (x, y, z) [m]
   */
  State AzElRangeToCart(const State& aer, const State& xyz_ref, Real R_body, Real flattening) {
    State enu = AzElRangeToEastNorthUp(aer);
    State xyz = EastNorthUpToCart(enu, xyz_ref, R_body, flattening);
    return Cart3(xyz, aer.GetFrame());
  }
  VEC_IMP_VECTOR_VECTOR(AzElRangeToCart, 3);

  /**
   * @brief Convert latitude, longitude, and altitude to stereographic coordinates
   * @param lla Latitude, longitude, and altitude (lat, lon, alt) [deg, deg, m]
   * @param R_body Body radius [m]
   * @return Stereographic coordinates (x, y, alt) [m]
   */
  State LatLonAltToStereographic(const State& lla, Real R_body) {
    LatLonAlt lla_tmp({lla(0), lla(1), 0}, lla.GetFrame());
    Cart3 xyz_cart = LatLonAltToCart(lla_tmp, R_body, 0);
    Real x = xyz_cart(0), y = xyz_cart(1), z = xyz_cart(2);

    Real xs = 2 * R_body / (R_body - z) * x;
    Real ys = 2 * R_body / (R_body - z) * y;
    Real alt = lla(2);
    return Cart3({xs, ys, alt}, lla.GetFrame());
  }
  VEC_IMP_VECTOR_REAL(LatLonAltToStereographic, 3);

  /**
   * @brief Convert stereographic coordinates to latitude, longitude, and altitude
   * @param xya Stereographic coordinates (x, y, alt) [m]
   * @param R_body Body radius [m]
   * @return Latitude, longitude, and altitude [deg, deg, m]
   */
  State StereographicToLatLonAlt(const State& xya, Real R_body) {
    Real xs = xya(0), ys = xya(1), alt = xya(2);
    Real denom = xs * xs + ys * ys + 4 * R_body * R_body;
    Real num = xs * xs + ys * ys - 4 * R_body * R_body;
    Real x = 4 * R_body * R_body * xs / denom;
    Real y = 4 * R_body * R_body * ys / denom;
    Real z = (num * R_body) / denom;

    LatLonAlt lla = CartToLatLonAlt(Cart3({x, y, z}, xya.GetFrame()), R_body, 0);
    lla(2) = alt;
    return lla;
  }
  VEC_IMP_VECTOR_REAL(StereographicToLatLonAlt, 3);

  /**
   * @brief Convert stereographic coordinates to Cartesian coordinates
   * @param xya Stereographic coordinates (x, y, alt) [m]
   * @param R_body Body radius [m]
   * @return Cartesian coordinates (x, y, z) [m]
   */
  State StereographicToCart(const State& xya, Real R_body) {
    LatLonAlt lla = StereographicToLatLonAlt(xya, R_body);
    return LatLonAltToCart(lla, R_body, 0);
  }
  VEC_IMP_VECTOR_REAL(StereographicToCart, 3);

  /**
   * @brief Convert Cartesian coordinates to stereographic coordinates
   * @param xyz Cartesian coordinates (x, y, z) [m]
   * @param R_body Body radius [m]
   * @return Stereographic coordinates (x, y, alt) [m]
   */
  State CartToStereographic(const State& xyz, Real R_body) {
    LatLonAlt lla = CartToLatLonAlt(xyz, R_body, 0);
    return LatLonAltToStereographic(lla, R_body);
  }
  VEC_IMP_VECTOR_REAL(CartToStereographic, 3);

}  // namespace lupnt
