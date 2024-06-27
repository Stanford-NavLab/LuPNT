#include "coordinates.h"

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
namespace lupnt {
Vec3 Spherical2Cart(const Vec3 &r_sph) {
  auto [r, theta, phi] = unpack(r_sph);
  Real x = r * sin(theta) * cos(phi);
  Real y = r * sin(theta) * sin(phi);
  Real z = r * cos(theta);
  return Vec3(x, y, z);
}

Vec3 Cart2Spherical(const Vec3 &r_cart) {
  auto [x, y, z] = unpack(r_cart);
  Real r = r_cart.norm();
  Real theta = acos(z / r);
  Real phi = atan2(y, x);
  return Vec3(r, theta, phi);
}

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

Vec3 Cart2LatLonAlt(const Vec3 &cart, Real R_body, Real flattening) {
  if (R_body == 0) R_body = cart.norm();

  auto [x, y, z] = unpack(cart);
  Real e2 = flattening * (2 - flattening);
  int max_iterations = 1000;

  Real delta_z = e2 * z;  // Initial guess for Δz
  Real delta_z_ant = delta_z + 1e3 * EPS;
  Real sin_phi, N;
  int it = 0;
  while (it < max_iterations && abs(delta_z_ant - delta_z) >= EPS) {
    delta_z_ant = delta_z;
    sin_phi =
        (z + delta_z) / sqrt(x * x + y * y + (z + delta_z) * (z + delta_z));
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

/// @brief
/// @param enu
/// @param r_ref
/// @param R_body
/// @param flattening
/// @return
Vec3 EastNorthUp2Cart(const Vec3 &enu, const Vec3 &r_ref, Real R_body,
                      Real flattening) {
  Vec3 lla = Cart2LatLonAlt(r_ref, R_body, flattening);
  auto [lat, lon, alt] = unpack(lla);
  Mat3 R = RotY(-lat) * RotZ(lon);
  Vec3 uen(enu(2), enu(0), enu(1));
  return r_ref + R.transpose() * uen;
}

/// @brief
/// @param r Cartesian position vector
/// @param r_ref Reference position vector
/// @param R_body Body radius
/// @param flattening Body flattening
/// @return
Vec3 Cart2EastNorthUp(const Vec3 &r, const Vec3 &r_ref, Real R_body,
                      Real flattening) {
  Vec3 lla = Cart2LatLonAlt(r_ref, R_body, flattening);
  Real lat = lla(0), lon = lla(1);
  Mat3 R = RotY(-lat) * RotZ(lon);
  Vec3 uen = R * (r - r_ref);
  Vec3 enu(uen(1), uen(2), uen(0));
  return enu;
}

Vec3 Cart2AzElRange(const Vec3 &r, const Vec3 &r_ref, Real R_body,
                    Real flattening) {
  Vec3 enu = Cart2EastNorthUp(r, r_ref, R_body, flattening);
  Real range = enu.norm();
  Real az = Wrap2TwoPi(atan2(enu(0), enu(1)));
  Real el = asin(enu(2) / range);
  return Vec3(az, el, range);
}

Vec3 AzElRange2Cart(const Vec3 &aer, const Vec3 &r_ref, Real R_body,
                    Real flattening) {
  auto [az, el, range] = unpack(aer);
  Real x = range * cos(el) * sin(az);
  Real y = range * cos(el) * cos(az);
  Real z = range * sin(el);
  Vec3 enu(x, y, z);
  return EastNorthUp2Cart(enu, r_ref, R_body, flattening);
}

VEC_IMP_VECTOR(Spherical2Cart, 3);
VEC_IMP_VECTOR(Cart2Spherical, 3);

}  // namespace lupnt