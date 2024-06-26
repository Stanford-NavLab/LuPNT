#include "axes_transformation.h"

namespace lupnt {
// Rotation Angles to DCM
Matrix3 angle2dcm(real roll, real pitch, real yaw) {
  Matrix3 dcm;
  dcm = AngleAxis(yaw, Vector3::UnitZ()) * AngleAxis(pitch, Vector3::UnitY()) *
        AngleAxis(roll, Vector3::UnitX());
  return dcm;
}

// Rotation Angles to Quaternion
Quaternion angle2quat(real roll, real pitch, real yaw) {
  Quaternion quat;
  quat = AngleAxis(yaw, Vector3::UnitZ()) * AngleAxis(pitch, Vector3::UnitY()) *
         AngleAxis(roll, Vector3::UnitX());
  return quat;
}

// Rotation Angles to Euler-Rodrigues Vector
Vector4 angle2rod(real roll, real pitch, real yaw) {
  Quaternion quat = angle2quat(roll, pitch, yaw);
  Vector4 rod;
  rod << quat.vec() * quat.w(), quat.vec();
  return rod;
}

// Euler-Rodrigues Vector to Rotation Angles
std::tuple<real, real, real> rod2angle(const Vector4 &rod) {
  Quaternion quat;
  quat.w() = sqrt(1 - rod.head<3>().squaredNorm());
  quat.vec() = rod.head<3>();
  return quat2angle(quat);
}

// Quaternion to Rotation Angles
std::tuple<real, real, real> quat2angle(const Quaternion &quat) {
  real roll, pitch, yaw;
  roll = atan2(2 * (quat.w() * quat.x() + quat.y() * quat.z()),
               1 - 2 * (quat.x() * quat.x() + quat.y() * quat.y()));
  pitch = asin(2 * (quat.w() * quat.y() - quat.z() * quat.x()));
  yaw = atan2(2 * (quat.w() * quat.z() + quat.x() * quat.y()),
              1 - 2 * (quat.y() * quat.y() + quat.z() * quat.z()));
  return std::make_tuple(roll, pitch, yaw);
}

// DCM to Angle of Attack and Sideslip Angle
std::tuple<real, real> dcm2alphabeta(const Matrix3 &dcm) {
  real alpha = atan2(dcm(2, 0), dcm(0, 0));
  real beta = asin(dcm(1, 0));
  return std::make_tuple(alpha, beta);
}

// DCM to Rotation Angles
std::tuple<real, real, real> dcm2angle(const Matrix3 &dcm) {
  real roll = atan2(dcm(2, 1), dcm(2, 2));
  real pitch = asin(-dcm(2, 0));
  real yaw = atan2(dcm(1, 0), dcm(0, 0));
  return std::make_tuple(roll, pitch, yaw);
}

// DCM to Geodetic Latitude and Longitude
std::tuple<real, real> dcm2latlon(const Matrix3 &dcm) {
  real lat = asin(dcm(2, 0));
  real lon = atan2(dcm(2, 1), dcm(2, 2));
  return std::make_tuple(lat, lon);
}

// Convert Body Frame to Stability Frame Transformation Matrix
Matrix3 dcmbody2stability(const Matrix3 &dcm) {
  // Assuming stability frame transformation is identity for this example
  return dcm;
}

// DCM to Quaternion
Quaternion dcm2quat(const Matrix3 &dcm) {
  Quaternion quat(dcm);
  return quat;
}

// DCM to Euler-Rodrigues Vector
Vector4 dcm2rod(const Matrix3 &dcm) {
  Quaternion quat = dcm2quat(dcm);
  Vector4 rod;
  rod << quat.vec() * quat.w(), quat.vec();
  return rod;
}

// Convert Angle of Attack and Sideslip Angle to DCM
Matrix3 dcmbody2wind(real angle_of_attack, real sideslip_angle) {
  Matrix3 dcm;
  dcm = AngleAxis(sideslip_angle, Vector3::UnitY()) *
        AngleAxis(angle_of_attack, Vector3::UnitX());
  return dcm;
}

// Convert Geodetic Latitude and Longitude to NED Direction Cosine Matrix
Matrix3 dcmecef2ned(real lat, real lon) {
  Matrix3 dcm;
  dcm = AngleAxis(-lon, Vector3::UnitZ()) * AngleAxis(-lat, Vector3::UnitY());
  return dcm;
}

// Convert ECI to ECEF Coordinates
Matrix3 dcmeci2ecef(real gmst) {
  Matrix3 dcm;
  dcm = AngleAxis(gmst, Vector3::UnitZ());
  return dcm;
}

// Quaternion to DCM
Matrix3 quat2dcm(const Quaternion &quat) { return quat.toRotationMatrix(); }

// Euler-Rodrigues Vector to DCM
Matrix3 rod2dcm(const Vector4 &rod) {
  Quaternion quat = rod2quat(rod);
  return quat.toRotationMatrix();
}

// Convert ECEF Coordinates to Geodetic Coordinates
Vector3 ecef2lla(const Vector3 &ecef) {
  // Placeholder implementation for ECEF to LLA conversion
  Vector3 lla;
  // Implement actual ECEF to LLA conversion logic here
  return lla;
}

// Convert ECI Coordinates to Latitude, Longitude, Altitude
Vector3 eci2lla(const Vector3 &eci, real gmst) {
  // Placeholder implementation for ECI to LLA conversion
  Vector3 lla;
  // Implement actual ECI to LLA conversion logic here
  return lla;
}

// Convert Flat Earth Position to Geodetic Coordinates
Vector3 flat2lla(const Vector3 &flat, real lat0, real lon0, real alt0) {
  // Placeholder implementation for Flat Earth to LLA conversion
  Vector3 lla;
  // Implement actual Flat Earth to LLA conversion logic here
  return lla;
}

// Convert Geodetic Coordinates to ECEF Coordinates
Vector3 lla2ecef(real lat, real lon, real alt) {
  real a = WGS84_A;
  real f = WGS84_F;
  real e2 = 2 * f - f * f;
  real sin_lat = sin(lat);
  real cos_lat = cos(lat);
  real N = a / sqrt(1 - e2 * sin_lat * sin_lat);
  real x = (N + alt) * cos_lat * cos(lon);
  real y = (N + alt) * cos_lat * sin(lon);
  real z = (N * (1 - e2) + alt) * sin_lat;
  return Vector3(x, y, z);
}

// Convert Geodetic Coordinates to ECI Coordinates
Vector3 lla2eci(real lat, real lon, real alt, real gmst) {
  Vector3 ecef = lla2ecef(lat, lon, alt);
  Matrix3 dcm = dcmeci2ecef(gmst).transpose();
  return dcm * ecef;
}

// Convert Geodetic Latitude, Longitude, and Altitude to Flat Earth Position
Vector3 lla2flat(real lat, real lon, real alt, real lat0, real lon0,
                 real alt0) {
  // Placeholder implementation for LLA to Flat Earth conversion
  Vector3 flat;
  // Implement actual LLA to Flat Earth conversion logic here
  return flat;
}

real geoc2geod(real geoc_lat) {
  real e2 = 2 * WGS84_F - WGS84_F * WGS84_F;
  return atan(tan(geoc_lat) / (1 - e2));
}

// Convert Geodetic Latitude to Geocentric Latitude
real geod2geoc(real geod_lat) {
  real e2 = 2 * WGS84_F - WGS84_F * WGS84_F;
  return atan((1 - e2) * tan(geod_lat));
}

// Convert Quaternion to Euler-Rodrigues Vector
Vector4 quat2rod(const Quaternion &quat) {
  Vector4 rod;
  rod << quat.vec() * quat.w(), quat.vec();
  return rod;
}

// Convert Quaternion to Rotation Angles
std::tuple<real, real, real> quat2angle(const Quaternion &quat) {
  real roll, pitch, yaw;
  roll = atan2(2 * (quat.w() * quat.x() + quat.y() * quat.z()),
               1 - 2 * (quat.x() * quat.x() + quat.y() * quat.y()));
  pitch = asin(2 * (quat.w() * quat.y() - quat.z() * quat.x()));
  yaw = atan2(2 * (quat.w() * quat.z() + quat.x() * quat.y()),
              1 - 2 * (quat.y() * quat.y() + quat.z() * quat.z()));
  return std::make_tuple(roll, pitch, yaw);
}

// Convert Euler-Rodrigues Vector to Quaternion
Quaternion rod2quat(const Vector4 &rod) {
  Quaternion quat;
  quat.w() = sqrt(1 - rod.head<3>().squaredNorm());
  quat.vec() = rod.head<3>();
  return quat;
}

// Convert Earth-centered inertial (ECI) coordinates to azimuth, elevation,
// slant range (AER) coordinates
std::tuple<real, real, real> eci2aer(const Vector3 &eci, real lat, real lon,
                                     real alt, real gmst) {
  Vector3 ecef = dcmeci2ecef(gmst) * eci;
  Vector3 lla = ecef2lla(ecef);

  real dlat = lla[0] - lat;
  real dlon = lla[1] - lon;
  real daz = atan2(sin(dlon) * cos(lla[0]),
                   cos(lat) * sin(lla[0]) - sin(lat) * cos(lla[0]) * cos(dlon));
  real del = asin(sin(lat) * sin(lla[0]) + cos(lat) * cos(lla[0]) * cos(dlon));
  real dsl = sqrt(eci.dot(eci));

  return std::make_tuple(daz, del, dsl);
}

// Convert direction cosine matrix to quaternion
Quaternion dcm2quat(const Matrix3 &dcm) {
  Quaternion quat(dcm);
  return quat;
}

// Convert direction cosine matrix to Euler-Rodrigues vector
Vector4 dcm2rod(const Matrix3 &dcm) {
  Quaternion quat = dcm2quat(dcm);
  Vector4 rod;
  rod << quat.vec() * quat.w(), quat.vec();
  return rod;
}

// Convert direction cosine matrix to rotation angles
std::tuple<real, real, real> dcm2angle(const Matrix3 &dcm) {
  real roll = atan2(dcm(2, 1), dcm(2, 2));
  real pitch = asin(-dcm(2, 0));
  real yaw = atan2(dcm(1, 0), dcm(0, 0));
  return std::make_tuple(roll, pitch, yaw);
}

// Convert direction cosine matrix to angle of attack and sideslip angle
std::tuple<real, real> dcm2alphabeta(const Matrix3 &dcm) {
  real alpha = atan2(dcm(2, 0), dcm(0, 0));
  real beta = asin(dcm(1, 0));
  return std::make_tuple(alpha, beta);
}

// Convert geodetic latitude and longitude to direction cosine matrix
Matrix3 dcmecef2ned(real lat, real lon) {
  Matrix3 dcm;
  dcm = AngleAxis(-lon, Vector3::UnitZ()) * AngleAxis(-lat, Vector3::UnitY());
  return dcm;
}

// Convert Earth-centered inertial (ECI) to Earth-centered Earth-fixed (ECEF)
// coordinates
Matrix3 dcmeci2ecef(real gmst) {
  Matrix3 dcm;
  dcm = AngleAxis(gmst, Vector3::UnitZ());
  return dcm;
}

// Convert from flat Earth position to array of geodetic coordinates
Vector3 flat2lla(const Vector3 &flat, real lat0, real lon0, real alt0) {
  real dlat = flat[0] / WGS84_A;
  real dlon = flat[1] / (WGS84_A * cos(lat0));
  real dalt = flat[2];
  return Vector3(lat0 + dlat, lon0 + dlon, alt0 + dalt);
}

// Convert from geodetic latitude, longitude, and altitude to flat Earth
// position
Vector3 lla2flat(real lat, real lon, real alt, real lat0, real lon0,
                 real alt0) {
  real dlat = (lat - lat0) * WGS84_A;
  real dlon = (lon - lon0) * WGS84_A * cos(lat0);
  real dalt = alt - alt0;
  return Vector3(dlat, dlon, dalt);
}

}  // namespace lupnt