/**
 * @file orbit_state_utils.h
 * @author Stanford NAV LAB
 * @brief Util functions for state conversions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "lupnt/core/constants.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/orbit_state.h"

// Function:
// Vec<size> = func(Vec<size>)
// New definitions:
// Mat<-1,size> = (Mat<-1,size>)
#define VEC_DEF_VECTOR(func, size) Mat<-1, size> func(const Mat<-1, size> &x);

// Function:
// Vec<size> = func(Vec<size>, real)
// New definitions:
// Mat<-1,size> = func(Vec<size>, VecX)
// Mat<-1,size> = func(Mat<-1,size>, real)
// Mat<-1,size> = func(Mat<-1,size>, VecX)
#define VEC_DEF_VECTOR_REAL(func, size)                  \
  Mat<-1, size> func(const Vec<size> &x, const VecX &y); \
  Mat<-1, size> func(const Mat<-1, size> &x, real y);    \
  Mat<-1, size> func(const Mat<-1, size> &x, const VecX &y);

// Function:
// Vec<size> = func(Vec<size>, real, real)
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, real, real)
#define VEC_DEF_VECTOR_REAL_REAL(func, size) \
  Mat<-1, size> func(const Mat<-1, size> &x, real y, real z);

// Function:
// Vec<size> = func(Vec<size>, Vec<size>, real)
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Vec<size>, real)
// Mat<-1,size> = func(Vec<size>, Mat<-1,size>, real)
// Mat<-1,size> = func(Mat<-1,size>, Mat<-1,size>, real)
#define VEC_DEF_VECTOR_VECTOR_REAL(func, size)                            \
  Mat<-1, size> func(const Mat<-1, size> &x, const Vec<size> &y, real z); \
  Mat<-1, size> func(const Vec<size> &x, const Mat<-1, size> &y, real z); \
  Mat<-1, size> func(const Mat<-1, size> &x, const Mat<-1, size> &y, real z);

// Function:
// Vec<size> = func(Vec<size>, Vec<size>
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Mat<-1,size>)
// Mat<-1,size> = func(Mat<-1,size>, Vec<size>)
// Mat<-1,size> = func(Vec<size>, Mat<-1,size>)
#define VEC_DEF_VECTOR_VECTOR(func, size)                             \
  Mat<-1, size> func(const Mat<-1, size> &x, const Mat<-1, size> &y); \
  Mat<-1, size> func(const Mat<-1, size> &x, const Vec<size> &y);     \
  Mat<-1, size> func(const Vec<size> &x, const Mat<-1, size> &y);

// Function:
// real = func(real, real)
// New definitions:
// vector = func(vector, real)
// vector = func(real, vector)
// vector = func(vector, vector)
#define VEC_DEF_REAL_REAL(func)     \
  VecX func(const VecX &x, real y); \
  VecX func(real x, const VecX &y); \
  VecX func(const VecX &x, const VecX &y);

namespace lupnt {

extern std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
                std::function<Vec6(const Vec6 &, real)>>
    absolute_conversions;

extern std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
                std::function<Vec6(const Vec6 &, const Vec6 &)>>
    relative_conversions;

Vec6 ConvertOrbitState(const Vec6 &state_in, OrbitStateRepres repres_in,
                       OrbitStateRepres repres_out, real mu);

Vec6 ConvertOrbitState(const Vec6 &state_in_c, const Vec6 &state_in_d,
                       OrbitStateRepres repres_in_c,
                       OrbitStateRepres repres_in_d,
                       OrbitStateRepres repres_out, real mu);

std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
    const std::shared_ptr<OrbitState> &state_in, OrbitStateRepres repres_out,
    real mu);

// From CartesianOrbitState
// - To ClassicalOE
ClassicalOE CartesianToClassical(const CartesianOrbitState &rv, real mu);
Vec6 CartesianToClassical(const Vec6 &rv, real mu);

// - To CartesianOrbitState (relative)
CartesianOrbitState InertialToRtn(const CartesianOrbitState &rv_c,
                                  const CartesianOrbitState &rv_d);
Vec6 InertialToRtn(const Vec6 &rv_c, const Vec6 &rv_d);

CartesianOrbitState RtnToInertial(const CartesianOrbitState &rv_c,
                                  const CartesianOrbitState &rv_rtn_d);
Vec6 RtnToInertial(const Vec6 &rv_c, const Vec6 &rv_rtn_d);

// From ClassicalOE
// - To CartesianOrbitState
CartesianOrbitState ClassicalToCartesian(const ClassicalOE &coe, real mu);
Vec6 ClassicalToCartesian(const Vec6 &coe, real mu);

// - To QuasiNonsingularOE
QuasiNonsingularOE ClassicalToQuasiNonsingular(const ClassicalOE &coe, real mu);
Vec6 ClassicalToQuasiNonsingular(const Vec6 &coe, real mu);

// - To EquinoctialOE
EquinoctialOE ClassicalToEquinoctial(const ClassicalOE &coe, real mu);
Vec6 ClassicalToEquinoctial(const Vec6 &coe, real mu);

// - To DelaunayOE
DelaunayOE ClassicalToDelaunay(const ClassicalOE &coe, real mu);
Vec6 ClassicalToDelaunay(const Vec6 &coe, real mu);

// From QuasiNonsingularOE
// - To ClassicalOE
ClassicalOE QuasiNonsingularToClassical(const QuasiNonsingularOE &qnsoe,
                                        real mu);
Vec6 QuasiNonsingularToClassical(const Vec6 &qnsoeVec, real mu);

// From EquinoctialOE
// - To ClassicalOE
ClassicalOE EquinoctialToClassical(const EquinoctialOE &eqoe, real mu);
Vec6 EquinoctialToClassical(const Vec6 &eqoe, real mu);

// From DelaunayOE
// - To ClassicalOE
ClassicalOE DelaunayToClassical(const DelaunayOE &deloe, real mu);
Vec6 DelaunayToClassical(const Vec6 &deloe, real mu);

// From ClassicalOE and QuasiNonsingularROE (relative)
// - To ClassicalOE
ClassicalOE RelativeQuasiNonsingularToClassical(
    const ClassicalOE &coe,
    const QuasiNonsingularROE &RelativeQuasiNonsingular);
Vec6 RelativeQuasiNonsingularToClassical(const Vec6 &coe,
                                         const Vec6 &RelativeQuasiNonsingular);

// Mean and Osculating
ClassicalOE OsculatingToMean(const ClassicalOE &coe_o, real J2);
ClassicalOE MeanToOsculating(const ClassicalOE &coe_m, real J2);
Vec6 OsculatingToMean(const Vec6 &coe_o, real mu, real J2);
Vec6 MeanToOsculating(const Vec6 &coe_m, real mu, real J2);

std::array<double, 6> ComputeSecondOrderShortPeriod(Vec6 &coe, Vec6 &doe);
std::array<double, 6> ComputeFirstOrderMediumPeriod(Vec6 &coe, Vec6 &doe);
std::array<double, 6> ComputeSecondOrderMediumPeriod(Vec6 &coe, Vec6 &doe);
std::array<double, 6> ComputeCorrectionMediumPeriod(Vec6 &coe, Vec6 &doe);

// Anomaly
real EccentricToTrueAnomaly(real E, real e);
real EccentricToMeanAnomaly(real E, real e);
real MeanToEccentricAnomaly(real M, real e);
real TrueToEccentricAnomaly(real nu, real e);
real MeanToTrueAnomaly(real M, real e);
real TrueToMeanAnomaly(real f, real e);

// Other coordinates
Vec3 LatLonAltToEcef(const Vec3 &r_geo, real radius, real flattening = 0);
Vec3 EcefToLatLonAlt(const Vec3 &r_cart, real radius, real flattening = 0);
Vec3 SphericalToCartesian(const Vec3 &r_sph);
Vec3 CartesianToSpherical(const Vec3 &r_cart);
Vec3 EastNorthUpToCartesian(const Vec3 &r_ref, const Vec3 &r_enu,
                            real flattening = -1.0);
Vec3 CartesianToEastNorthUp(const Vec3 &r_ref, const Vec3 &r_cart,
                            real flattening = -1.0);
Vec3 CartesianToAzimuthElevationRange(const Vec3 &r_cart_ref,
                                      const Vec3 &r_cart,
                                      real flattening = -1.0);
Vec3 AzimuthElevationRangeToCartesian(const Vec3 &r_aer_ref, const Vec3 &r_aer,
                                      real flattening = -1.0);

VEC_DEF_REAL_REAL(EccentricToTrueAnomaly);
VEC_DEF_REAL_REAL(EccentricToMeanAnomaly);
VEC_DEF_REAL_REAL(MeanToEccentricAnomaly);
VEC_DEF_REAL_REAL(TrueToEccentricAnomaly);
VEC_DEF_REAL_REAL(MeanToTrueAnomaly);
VEC_DEF_REAL_REAL(TrueToMeanAnomaly);

VEC_DEF_VECTOR(SphericalToCartesian, 3);
VEC_DEF_VECTOR(CartesianToSpherical, 3);
VEC_DEF_VECTOR_REAL(LatLonAltToEcef, 3);
VEC_DEF_VECTOR_REAL(EcefToLatLonAlt, 3);
VEC_DEF_VECTOR_REAL_REAL(GeodeticToCartesian, 3);
VEC_DEF_VECTOR_REAL_REAL(CartesianToGeodetic, 3);

VEC_DEF_VECTOR_VECTOR(EastNorthUpToCartesian, 3);
VEC_DEF_VECTOR_VECTOR(CartesianToEastNorthUp, 3);
VEC_DEF_VECTOR_VECTOR(CartesianToAzimuthElevationRange, 3);
VEC_DEF_VECTOR_VECTOR(AzimuthElevationRangeToCartesian, 3);
VEC_DEF_VECTOR_VECTOR_REAL(EastNorthUpToCartesian, 3);
VEC_DEF_VECTOR_VECTOR_REAL(CartesianToEastNorthUp, 3);
VEC_DEF_VECTOR_VECTOR_REAL(CartesianToAzimuthElevationRange, 3);
VEC_DEF_VECTOR_VECTOR_REAL(AzimuthElevationRangeToCartesian, 3);

VEC_DEF_VECTOR_REAL(ClassicalToCartesian, 6);
VEC_DEF_VECTOR_VECTOR(InertialToRtn, 6);
VEC_DEF_VECTOR_VECTOR(RtnToInertial, 6);
VEC_DEF_VECTOR_REAL(CartesianToClassical, 6);
VEC_DEF_VECTOR_REAL(ClassicalToQuasiNonsingular, 6);
VEC_DEF_VECTOR_REAL(ClassicalToEquinoctial, 6);
VEC_DEF_VECTOR_REAL(ClassicalToDelaunay, 6);
VEC_DEF_VECTOR_REAL(QuasiNonsingularToClassical, 6);
VEC_DEF_VECTOR_REAL(EquinoctialToClassical, 6);
VEC_DEF_VECTOR_REAL(DelaunayToClassical, 6);
VEC_DEF_VECTOR_VECTOR(RelativeQuasiNonsingularToClassical, 6);

class TLE {
 public:
  std::string name;
  double epochYear;
  double epochDay;
  double epochTAI;
  double bstar;
  double inclination;
  double raan;
  double eccentricity;
  double argPerigee;
  double meanAnomaly;
  double meanMotion;
  int prn;

  static TLE FromLines(const std::string &line1, const std::string &line2,
                       const std::string &line3);
  static std::vector<TLE> FromFile(const std::string &filename);
};

static std::shared_ptr<CartesianOrbitState> ConvertOrbitStateFrame(
    const std::shared_ptr<CartesianOrbitState> state_in, const real epoch,
    const Frame frame_out) {
  auto rv_in = state_in->GetVec();
  auto frame_in = state_in->GetCoordSystem();
  auto rv_out = FrameConverter::Convert(epoch, rv_in, frame_in, frame_out);
  auto state_out = std::make_shared<CartesianOrbitState>(rv_out, frame_out);
  return state_out;
}

}  // namespace lupnt