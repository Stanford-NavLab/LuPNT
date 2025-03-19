#include "lupnt/physics/coordinates.h"

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
namespace lupnt {

  /// @brief Convert spherical coordinates to Cartesian coordinates
  /// @param enu East-North-Up coordinates
  /// @return Cartesian coordinates
  Vec3 EastNorthUp2AzElRange(const Vec3 &enu) {
    auto [e, n, u] = unpack(enu);
    Real range = enu.norm();
    Real azim = atan2(e, n);
    Real elev = asin(u / range);
    return Vec3(azim, elev, range);
  }

  /// @brief Convert Cartesian coordinates to spherical coordinates
  /// @param aer Azimuth-Elevation-Range coordinates
  /// @return Cartesian coordinates
  Vec3 AzElRange2EastNorthUp(const Vec3 &aer) {
    auto [az, el, range] = unpack(aer);
    Real x = range * cos(el) * sin(az);
    Real y = range * cos(el) * cos(az);
    Real z = range * sin(el);
    Vec3 enu(x, y, z);
    return enu;
  }

  /// @brief Convert latitude, longitude, and altitude to Cartesian coordinates
  /// @param lla Latitude, longitude, and altitude
  /// @param R_body Body radius
  /// @param flattening Body flattening
  /// @return Cartesian coordinates
  Vec3 LatLonAlt2Cart(const Vec3 &lla, Real R_body, Real flattening) {
    auto [lat, lon, alt] = unpack(lla);
    if (R_body == 0) R_body = alt;

    Real e2 = flattening * (2 - flattening);
    Real cos_lat = cos(lat);
    Real sin_lat = sin(lat);
    Real N = R_body / sqrt(1 - e2 * sin_lat * sin_lat);
    Real x = (N + alt) * cos_lat * cos(lon);
    Real y = (N + alt) * cos_lat * sin(lon);
    Real z = ((1 - e2) * N + alt) * sin_lat;
    return Vec3(x, y, z);
  }

  /// @brief Convert Cartesian coordinates to latitude, longitude, and altitude
  /// @param cart Cartesian coordinates
  /// @param R_body Body radius
  /// @param flattening Body flattening
  /// @return Latitude, longitude, and altitude
  Vec3 Cart2LatLonAlt(const Vec3 &cart, Real R_body, Real flattening) {
    if (R_body == 0) R_body = cart.norm();

    auto [x, y, z] = unpack(cart);
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
    assert(it < max_iterations && "Cart2LatLonAlt did not converge");

    Real lon = atan2(y, x);
    Real lat = atan2(z + delta_z, sqrt(x * x + y * y));
    Real h = sqrt(x * x + y * y + (z + delta_z) * (z + delta_z)) - N;
    return Vec3(lat, lon, h);
  }

  /// @brief Convert East-North-Up coordinates to Cartesian coordinates
  /// @param enu East-North-Up coordinates
  /// @param xyz_ref Reference position vector
  /// @param R_body Body radius
  /// @param flattening Body flattening
  Vec3 EastNorthUp2Cart(const Vec3 &enu, const Vec3 &xyz_ref, Real R_body, Real flattening) {
    Vec3 lla = Cart2LatLonAlt(xyz_ref, R_body, flattening);
    auto [lat, lon, alt] = unpack(lla);
    Mat3 R = RotX(PI_OVER_TWO - lat) * RotZ(PI_OVER_TWO + lon);
    Vec3 xyz = R.transpose() * enu + xyz_ref;
    return xyz;
  }

  /// @brief
  /// @param xyz Cartesian position vector
  /// @param xyz_ref Reference position vector
  /// @param R_body Body radius
  /// @param flattening Body flattening
  /// @return
  Vec3 Cart2EastNorthUp(const Vec3 &xyz, const Vec3 &xyz_ref, Real R_body, Real flattening) {
    Vec3 lla = Cart2LatLonAlt(xyz_ref, R_body, flattening);
    auto [lat, lon, alt] = unpack(lla);
    Mat3 R = RotX(PI_OVER_TWO - lat) * RotZ(PI_OVER_TWO + lon);
    Vec3 enu = R * (xyz - xyz_ref);
    return enu;
  }

  /// @brief Convert Cartesian coordinates to Azimuth-Elevation-Range coordinates
  /// @param xyz Cartesian position vector
  /// @param xyz_ref Reference position vector
  /// @param R_body Body radius
  /// @param flattening Body flattening
  Vec3 Cart2AzElRange(const Vec3 &xyz, const Vec3 &xyz_ref, Real R_body, Real flattening) {
    Vec3 enu = Cart2EastNorthUp(xyz, xyz_ref, R_body, flattening);
    Vec3 aer = EastNorthUp2AzElRange(enu);
    return aer;
  }

  /// @brief Convert Azimuth-Elevation-Range coordinates to Cartesian coordinates
  /// @param aer Azimuth-Elevation-Range coordinates
  /// @param xyz_ref Reference position vector
  /// @param R_body Body radius
  /// @param flattening Body flattening
  Vec3 AzElRange2Cart(const Vec3 &aer, const Vec3 &xyz_ref, Real R_body, Real flattening) {
    Vec3 enu = AzElRange2EastNorthUp(aer);
    Vec3 xyz = EastNorthUp2Cart(enu, xyz_ref, R_body, flattening);
    return xyz;
  }

  // Vector implementations
  VEC_IMP_VECTOR(LatLonAlt2Cart, 3)
  VEC_IMP_VECTOR(Cart2LatLonAlt, 3)
  VEC_IMP_VECTOR_REAL(LatLonAlt2Cart, 3)
  VEC_IMP_VECTOR_REAL(Cart2LatLonAlt, 3)
  VEC_IMP_VECTOR_REAL_REAL(LatLonAlt2Cart, 3)
  VEC_IMP_VECTOR_REAL_REAL(Cart2LatLonAlt, 3)

  VEC_IMP_VECTOR_VECTOR(EastNorthUp2Cart, 3)
  VEC_IMP_VECTOR_VECTOR(Cart2EastNorthUp, 3)

  VEC_IMP_VECTOR_VECTOR(Cart2AzElRange, 3)
  VEC_IMP_VECTOR_VECTOR(AzElRange2Cart, 3)
}  // namespace lupnt
