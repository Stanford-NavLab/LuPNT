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
#include "lupnt/physics/coord_converter.h"
#include "lupnt/physics/orbit_state.h"

// Function:
// Vector<size> = func(Vector<size>)
// New definitions:
// Matrix<-1,size> = (Matrix<-1,size>)
#define VECTORIZED_DEFINITION_FROM_VECTOR(func, size) \
  Matrix<-1, size> func(const Matrix<-1, size> &x);

// Function:
// Vector<size> = func(Vector<size>, real)
// New definitions:
// Matrix<-1,size> = func(Matrix<-1,size>, real)
// Matrix<-1,size> = func(Vector<size>, VectorX)
// Matrix<-1,size> = func(Matrix<-1,size>, VectorX)
#define VECTORIZED_DEFINITION_FROM_VECTOR_REAL(func, size)        \
  Matrix<-1, size> func(const Vector<size> &x, const VectorX &y); \
  Matrix<-1, size> func(const Matrix<-1, size> &x, real y);       \
  Matrix<-1, size> func(const Matrix<-1, size> &x, const VectorX &y);

// Function:
// Vector<size> = func(Vector<size>, Vector<size>
// New definitions:
// Matrix<-1,size> = func(Matrix<-1,size>, Matrix<-1,size>)
// Matrix<-1,size> = func(Matrix<-1,size>, Vector<size>)
// Matrix<-1,size> = func(Vector<size>, Matrix<-1,size>)
#define VECTORIZED_DEFINITION_FROM_VECTOR_VECTOR(func, size)                   \
  Matrix<-1, size> func(const Matrix<-1, size> &x, const Matrix<-1, size> &y); \
  Matrix<-1, size> func(const Matrix<-1, size> &x, const Vector<size> &y);     \
  Matrix<-1, size> func(const Vector<size> &x, const Matrix<-1, size> &y);

// Function:
// real = func(real, real)
// New definitions:
// vector = func(vector, real)
// vector = func(real, vector)
// vector = func(vector, vector)
#define VECTORIZED_DEFINITION_FROM_REAL_REAL(func) \
  VectorX func(const VectorX &x, real y);          \
  VectorX func(real x, const VectorX &y);          \
  VectorX func(const VectorX &x, const VectorX &y);

namespace lupnt {

extern std::map<std::tuple<OrbitStateRepres, OrbitStateRepres>,
                std::function<Vector6(const Vector6 &, double)>>
    absolute_conversions;

extern std::map<std::tuple<OrbitStateRepres, OrbitStateRepres>,
                std::function<Vector6(const Vector6 &, const Vector6 &)>>
    relative_conversions;

Vector6 ConvertOrbitState(const Vector6 &state_in, OrbitStateRepres repres_in,
                          OrbitStateRepres repres_out, double mu = MU_MOON);

Vector6 ConvertOrbitState(const Vector6 &state_in_c, const Vector6 &state_in_d,
                          OrbitStateRepres repres_in_c,
                          OrbitStateRepres repres_in_d,
                          OrbitStateRepres repres_out, double mu = MU_MOON);

std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
    const std::shared_ptr<OrbitState> &state_in, OrbitStateRepres repres_out,
    double mu = MU_MOON);

// From CartesianOrbitState
// - To ClassicalOE
ClassicalOE CartesianToClassical(const CartesianOrbitState &rv,
                                 double mu = MU_MOON);
Vector6 CartesianToClassical(const Vector6 &rv, double mu = MU_MOON);

// - To CartesianOrbitState (relative)
CartesianOrbitState InertialToRtn(const CartesianOrbitState &rv_c,
                                  const CartesianOrbitState &rv_d);
Vector6 InertialToRtn(const Vector6 &rv_c, const Vector6 &rv_d);

CartesianOrbitState RtnToInertial(const CartesianOrbitState &rv_c,
                                  const CartesianOrbitState &rv_rtn_d);
Vector6 RtnToInertial(const Vector6 &rv_c, const Vector6 &rv_rtn_d);

// From ClassicalOE
// - To CartesianOrbitState
CartesianOrbitState ClassicalToCartesian(const ClassicalOE &coe,
                                         double mu = MU_MOON);
Vector6 ClassicalToCartesian(const Vector6 &coe, double mu = MU_MOON);

// - To QuasiNonsingularOE
QuasiNonsingularOE ClassicalToQuasiNonsingular(const ClassicalOE &coe,
                                               double mu = MU_MOON);
Vector6 ClassicalToQuasiNonsingular(const Vector6 &coe, double mu = MU_MOON);

// - To EquinoctialOE
EquinoctialOE ClassicalToEquinoctial(const ClassicalOE &coe,
                                     double mu = MU_MOON);
Vector6 ClassicalToEquinoctial(const Vector6 &coe, double mu = MU_MOON);

// - To DelaunayOE
DelaunayOE ClassicalToDelaunay(const ClassicalOE &coe, double mu = MU_MOON);
Vector6 ClassicalToDelaunay(const Vector6 &coe, double mu = MU_MOON);

// From QuasiNonsingularOE
// - To ClassicalOE
ClassicalOE QuasiNonsingularToClassical(const QuasiNonsingularOE &qnsoe,
                                        double mu = MU_MOON);
Vector6 QuasiNonsingularToClassical(const Vector6 &qnsoeVec,
                                    double mu = MU_MOON);

// From EquinoctialOE
// - To ClassicalOE
ClassicalOE EquinoctialToClassical(const EquinoctialOE &eqoe,
                                   double mu = MU_MOON);
Vector6 EquinoctialToClassical(const Vector6 &eqoe, double mu = MU_MOON);

// From DelaunayOE
// - To ClassicalOE
ClassicalOE DelaunayToClassical(const DelaunayOE &deloe, double mu = MU_MOON);
Vector6 DelaunayToClassical(const Vector6 &deloe, double mu = MU_MOON);

// From ClassicalOE and QuasiNonsingularROE (relative)
// - To ClassicalOE
ClassicalOE RelativeQuasiNonsingularToClassical(
    const ClassicalOE &coe,
    const QuasiNonsingularROE &RelativeQuasiNonsingular);
Vector6 RelativeQuasiNonsingularToClassical(
    const Vector6 &coe, const Vector6 &RelativeQuasiNonsingular);

// Mean and Osculating
ClassicalOE OsculatingToMean(const ClassicalOE &coe_o, double J2);
ClassicalOE MeanToOsculating(const ClassicalOE &coe_m, double J2);
Vector6 OsculatingToMean(const Vector6 &coe_o, double J2);
Vector6 MeanToOsculating(const Vector6 &coe_m, double J2);

std::array<double, 6> ComputeSecondOrderShortPeriod(Vector6 &coe, Vector6 &doe);
std::array<double, 6> ComputeFirstOrderMediumPeriod(Vector6 &coe, Vector6 &doe);
std::array<double, 6> ComputeSecondOrderMediumPeriod(Vector6 &coe,
                                                     Vector6 &doe);
std::array<double, 6> ComputeCorrectionMediumPeriod(Vector6 &coe, Vector6 &doe);

// Anomaly
real EccentricToTrueAnomaly(real E, real e);
real EccentricToMeanAnomaly(real E, real e);
real MeanToEccentricAnomaly(real M, real e);
real TrueToEccentricAnomaly(real nu, real e);
real MeanToTrueAnomaly(real M, real e);
real TrueToMeanAnomaly(real f, real e);

// Other coordinates
Vector3 GeographicalToCartesian(const Vector3 &r_geo, real radius);
Vector3 CartesianToGeographical(const Vector3 &r_cart, real radius);
Vector3 SphericalToCartesian(const Vector3 &r_sph);
Vector3 CartesianToSpherical(const Vector3 &r_cart);
Vector3 EastNorthUpToCartesian(const Vector3 &r_ref, const Vector3 &r_enu);
Vector3 CartesianToEastNorthUp(const Vector3 &r_ref, const Vector3 &r_cart);
Vector3 CartesianToAzimuthElevationRange(const Vector3 &r_cart_ref,
                                         const Vector3 &r_cart);
Vector3 AzimuthElevationRangeToCartesian(const Vector3 &r_aer_ref,
                                         const Vector3 &r_aer);

VECTORIZED_DEFINITION_FROM_REAL_REAL(EccentricToTrueAnomaly);
VECTORIZED_DEFINITION_FROM_REAL_REAL(EccentricToMeanAnomaly);
VECTORIZED_DEFINITION_FROM_REAL_REAL(MeanToEccentricAnomaly);
VECTORIZED_DEFINITION_FROM_REAL_REAL(TrueToEccentricAnomaly);
VECTORIZED_DEFINITION_FROM_REAL_REAL(MeanToTrueAnomaly);
VECTORIZED_DEFINITION_FROM_REAL_REAL(TrueToMeanAnomaly);

VECTORIZED_DEFINITION_FROM_VECTOR(SphericalToCartesian, 3);
VECTORIZED_DEFINITION_FROM_VECTOR(CartesianToSpherical, 3);
VECTORIZED_DEFINITION_FROM_VECTOR_REAL(GeographicalToCartesian, 3);
VECTORIZED_DEFINITION_FROM_VECTOR_REAL(CartesianToGeographical, 3);
VECTORIZED_DEFINITION_FROM_VECTOR_VECTOR(EastNorthUpToCartesian, 3);
VECTORIZED_DEFINITION_FROM_VECTOR_VECTOR(CartesianToEastNorthUp, 3);
VECTORIZED_DEFINITION_FROM_VECTOR_VECTOR(CartesianToAzimuthElevationRange, 3);
VECTORIZED_DEFINITION_FROM_VECTOR_VECTOR(AzimuthElevationRangeToCartesian, 3);

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

static std::shared_ptr<CartesianOrbitState> ConvertOrbitStateCoordSystem(
    const std::shared_ptr<CartesianOrbitState> state_in, const real epoch,
    const CoordSystem coord_sys_out) {
  auto rv_in = state_in->GetVector();
  auto coord_sys_in = state_in->GetCoordSystem();
  auto rv_out =
      CoordConverter::Convert(epoch, rv_in, coord_sys_in, coord_sys_out);
  auto state_out = std::make_shared<CartesianOrbitState>(rv_out, coord_sys_out);
  return state_out;
}

}  // namespace lupnt