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
ClassicalOE CartToCoe(const CartesianOrbitState &cart, double mu);
Vector6real CoeToCart(const Vector6real &coe, double mu);
Vector6real CartToCoe(const Vector6real &cart, double mu);

// COE <-> ROE
ClassicalOE QnsroeToCoe(const ClassicalOE &coe,
                        const QuasiNonsingularROE &qnsroe);
Vector6real QnsroeToCoe(const Vector6real &coe, const Vector6real &qnsroe);

// Inertial <-> RTN
CartesianOrbitState InertialToRtn(const CartesianOrbitState &cart_c,
                                  const CartesianOrbitState &cart_d);
Vector6real InertialToRtn(const Vector6real &cart_c, const Vector6real &cart_d);

// COE <-> RTN
CartesianOrbitState CoeToRtn(const ClassicalOE &coe_c, const ClassicalOE &coe_d,
                             double mu);
Vector6real CoeToRtn(const Vector6real &coe_c, const Vector6real &coe_d,
                     double mu);

// COE <-> QNSOE
QuasiNonsingularOE CoeToQnsoe(const ClassicalOE &coe);
ClassicalOE QnsoeToCoe(const QuasiNonsingularOE &qnsoe);
Vector6real CoeToQnsoe(const Vector6real &coe);
Vector6real QnsoeToCoe(const Vector6real &qnsoeVec);

// COE <-> QNSROE
QuasiNonsingularROE QnsoeToQnsroe(const QuasiNonsingularOE &qnsoe_c,
                                  const QuasiNonsingularOE &qnsoe_d);
QuasiNonsingularROE CoeToQnsroe(const ClassicalOE &coe_c,
                                const ClassicalOE &coe_d);

// COE <-> Equinoctial
EquinoctialOE CoeToEqoe(const ClassicalOE &coe);
ClassicalOE EqoeToCoe(const EquinoctialOE &eqoe);
Vector6real CoeToEqoe(const Vector6real &coe);
Vector6real EqoeToCoe(const Vector6real &eqoe);

// Mean <-> Osculating
ClassicalOE OscToMean(const ClassicalOE &coe_o, double J2);
ClassicalOE MeanToOsc(const ClassicalOE &coe_m, double J2);
Vector6real OscToMean(const Vector6real &coe_o, double J2);
Vector6real MeanToOsc(const Vector6real &coe_m, double J2);

Vector6real CartToQnsoe(const Vector6real &cart, double mu);
QuasiNonsingularOE CartToQnsoe(const CartesianOrbitState &cart, double mu);

Vector6real DelaunayToCoe(const Vector6real &deloe, double mu, double n,
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

real EccentricToTrueAnomaly(real E, real e);
real EccentricToMeanAnomaly(real E, real e);
real MeanToEccentricAnomaly(real M, real e);
real TrueToEccentricAnomaly(real nu, real e);
real MeanToTrueAnomaly(real M, real e);
real TrueToMeanAnomaly(real f, real e);

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