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

#include "lupnt/core/constants.h"
#include "lupnt/physics/coord_converter.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

// COE <-> Cart
CartesianOrbitState CoeToCart(const ClassicalOE &coe, double mu);
Vector6real CoeToCart(const Vector6real &coeVec, double mu);

ClassicalOE CartToCoe(const CartesianOrbitState &cartOrbitState, double mu);
Vector6real CartToCoe(const Vector6real &cartVec, double mu);

// ROE <-> COE
ClassicalOE RoeToCoe(const ClassicalOE &coeChief,
                     const QuasiNonsingularROE &roe);
Vector6real RoeToCoe(const Vector6real &coeChiefVec, const Vector6real &roeVec);

// Inertial <-> RTN
CartesianOrbitState InertialToRtn(
    const CartesianOrbitState &rtnOrigin,
    const CartesianOrbitState &inertialOrbitState);
Vector6real InertialToRtn(const Vector6real &rtnOrigin,
                          const Vector6real &inertialVec);

// COE <-> RTN
CartesianOrbitState CoeToRtn(const ClassicalOE &coeChief,
                             const ClassicalOE &coeDeputy, double mu);
Vector6real CoeToRtn(const Vector6real &coeChiefVec,
                     const Vector6real &coeDeputyVec, double mu);

// COE <-> QNSOE
QuasiNonsingularOE CoeToQnsoe(const ClassicalOE &coe);
Vector6real CoeToQnsoe(const Vector6real &coeVec);

// COE <-> QNSROE
QuasiNonsingularROE QnsoeToQnsroe(const QuasiNonsingularOE &qnsoeChief,
                                  const QuasiNonsingularOE &qnsoeDeputy);
QuasiNonsingularROE CoeToQnsroe(const ClassicalOE &coeChief,
                                const ClassicalOE &coeDeputy);

ClassicalOE QnsoeToCoe(const QuasiNonsingularOE &qnsoe);
Vector6real QnsoeToCoe(const Vector6real &qnsoeVec);

Vector6real EquioeToCoe(const Vector6real &equioe);
Vector6real CoeToEquioe(const Vector6real &coe);

Vector6real MeanToOsculating(const Vector6real &meanCoe, double J2);
ClassicalOE MeanToOsculating(const ClassicalOE &meanCoe, double J2);

Vector6real OsculatingToMean(const Vector6real &osculatingCoe, double J2);
ClassicalOE OsculatingToMean(const ClassicalOE &osculatingCoe, double J2);

Vector6real CartToQnsoe(const Vector6real &cart, double mu);
QuasiNonsingularOE CartToQnsoe(const CartesianOrbitState &cart, double mu);

QuasiNonsingularOE OscQnsoeToMeanQnsoe(const QuasiNonsingularOE &oscQnsoe,
                                       double J2);

Vector6real DelaunayToCoe(const Vector6real &delaunay, double mu, double n,
                          double t);
Vector6real CoeToDelaunay(const Vector6real &coe, double mu, double n,
                          double t);

std::array<double, 6> ComputeSecondOrderShortPeriod(Vector6real &coe,
                                                    Vector6real &doe);
std::array<double, 6> ComputeFirstOrderMediumPeriod(Vector6real &coe,
                                                    Vector6real &doe);
std::array<double, 6> ComputeSecondOrderMediumPeriod(Vector6real &coe,
                                                     Vector6real &doe);
std::array<double, 6> ComputeCorrectionMediumPeriod(Vector6real &coe,
                                                    Vector6real &doe);

// Anomaly conversions

real EccentricAnomToTrueAnom(real E, real e);
real EccentricAnomToMeanAnom(real E, real e);
real MeanAnomToEccentricAnom(real M, real e);
real MeanAnomToTrueAnom(real M, real e);
real TrueAnomToEccentricAnom(real nu, real e);
real TrueAnomToMeanAnom(real f, real e);

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
      CoordConverter::Convert(rv_in, epoch, coord_sys_in, coord_sys_out);
  auto state_out = std::make_shared<CartesianOrbitState>(rv_out, coord_sys_out);
  return state_out;
}

}  // namespace lupnt