#pragma once
#include "lupnt/core/constants.h"

namespace lupnt {

// Rotation Angles
Mat3 Angles2Dcm(real ang1, real ang2, real ang3,
                const std::string &order = "ZYX");
Quat Angles2Quat(real ang1, real ang2, real ang3,
                 const std::string &order = "ZYX");
Vec4 Angles2Rodrigues(real ang1, real ang2, real ang3,
                      const std::string &order = "ZYX");
Vec3 Rodrigues2Angles(const Vec4 &rod, const std::string &order = "ZYX");
Vec3 Quat2Angles(const Quat &quat);

// Direction Cosine Mat
Mat3 angle2dcm(real roll, real pitch, real yaw);
std::tuple<real, real> dcm2alphabeta(const Mat3 &dcm);
std::tuple<real, real, real> dcm2angle(const Mat3 &dcm);
std::tuple<real, real> dcm2latlon(const Mat3 &dcm);
Mat3 dcmbody2stability(const Mat3 &dcm);
Quat dcm2quat(const Mat3 &dcm);
Vec4 dcm2rod(const Mat3 &dcm);
Mat3 dcmbody2wind(real angle_of_attack, real sideslip_angle);
Mat3 dcmecef2ned(real lat, real lon);
Mat3 dcmeci2ecef(real gmst);
Mat3 quat2dcm(const Quat &quat);
Mat3 rod2dcm(const Vec4 &rod);

// Geodetic Latitude, Longitude, and Altitude
Vec3 ecef2lla(const Vec3 &ecef);
Vec3 eci2lla(const Vec3 &eci, real gmst);
Vec3 flat2lla(const Vec3 &flat, real lat0, real lon0, real alt0);
Vec3 lla2ecef(real lat, real lon, real alt);
Vec3 lla2eci(real lat, real lon, real alt, real gmst);
Vec3 lla2flat(real lat, real lon, real alt, real lat0, real lon0, real alt0);

// Geocentric Latitude and Geodetic Latitude
real geoc2geod(real geoc_lat);
real geod2geoc(real geod_lat);

// Quaternions and Rotation Angles
Quat angle2quat(real roll, real pitch, real yaw);
Quat dcm2quat(const Mat3 &dcm);
std::tuple<real, real, real> quat2angle(const Quat &quat);
Mat3 quat2dcm(const Quat &quat);
Vec4 quat2rod(const Quat &quat);
Quat rod2quat(const Vec4 &rod);

// Earth-Centered Inertial (ECI) to Local Azimuth
std::tuple<real, real, real> eci2aer(const Vec3 &eci, real lat, real lon,
                                     real alt, real gmst);

// Rodrigues Vecs
Vec4 angle2rod(real roll, real pitch, real yaw);
Vec4 dcm2rod(const Mat3 &dcm);
Vec4 quat2rod(const Quat &quat);
std::tuple<real, real, real> rod2angle(const Vec4 &rod);
Mat3 rod2dcm(const Vec4 &rod);
Quat rod2quat(const Vec4 &rod);

}  // namespace lupnt