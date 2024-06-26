#pragma once
#include "lupnt/core/constants.h"

namespace lupnt {

// Rotation Angles
Matrix3 Angles2Dcm(real ang1, real ang2, real ang3,
                   const std::string &order = "ZYX");
Quaternion Angles2Quat(real ang1, real ang2, real ang3,
                       const std::string &order = "ZYX");
Vector4 Angles2Rodrigues(real ang1, real ang2, real ang3,
                         const std::string &order = "ZYX");
Vector3 Rodrigues2Angles(const Vector4 &rod, const std::string &order = "ZYX");
Vector3 Quat2Angles(const Quaternion &quat);

// Direction Cosine Matrix
Matrix3 angle2dcm(real roll, real pitch, real yaw);
std::tuple<real, real> dcm2alphabeta(const Matrix3 &dcm);
std::tuple<real, real, real> dcm2angle(const Matrix3 &dcm);
std::tuple<real, real> dcm2latlon(const Matrix3 &dcm);
Matrix3 dcmbody2stability(const Matrix3 &dcm);
Quaternion dcm2quat(const Matrix3 &dcm);
Vector4 dcm2rod(const Matrix3 &dcm);
Matrix3 dcmbody2wind(real angle_of_attack, real sideslip_angle);
Matrix3 dcmecef2ned(real lat, real lon);
Matrix3 dcmeci2ecef(real gmst);
Matrix3 quat2dcm(const Quaternion &quat);
Matrix3 rod2dcm(const Vector4 &rod);

// Geodetic Latitude, Longitude, and Altitude
Vector3 ecef2lla(const Vector3 &ecef);
Vector3 eci2lla(const Vector3 &eci, real gmst);
Vector3 flat2lla(const Vector3 &flat, real lat0, real lon0, real alt0);
Vector3 lla2ecef(real lat, real lon, real alt);
Vector3 lla2eci(real lat, real lon, real alt, real gmst);
Vector3 lla2flat(real lat, real lon, real alt, real lat0, real lon0, real alt0);

// Geocentric Latitude and Geodetic Latitude
real geoc2geod(real geoc_lat);
real geod2geoc(real geod_lat);

// Quaternions and Rotation Angles
Quaternion angle2quat(real roll, real pitch, real yaw);
Quaternion dcm2quat(const Matrix3 &dcm);
std::tuple<real, real, real> quat2angle(const Quaternion &quat);
Matrix3 quat2dcm(const Quaternion &quat);
Vector4 quat2rod(const Quaternion &quat);
Quaternion rod2quat(const Vector4 &rod);

// Earth-Centered Inertial (ECI) to Local Azimuth
std::tuple<real, real, real> eci2aer(const Vector3 &eci, real lat, real lon,
                                     real alt, real gmst);

// Rodrigues Vectors
Vector4 angle2rod(real roll, real pitch, real yaw);
Vector4 dcm2rod(const Matrix3 &dcm);
Vector4 quat2rod(const Quaternion &quat);
std::tuple<real, real, real> rod2angle(const Vector4 &rod);
Matrix3 rod2dcm(const Vector4 &rod);
Quaternion rod2quat(const Vector4 &rod);

}  // namespace lupnt