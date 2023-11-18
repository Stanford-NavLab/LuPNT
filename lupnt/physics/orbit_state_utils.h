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
#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <boost/preprocessor.hpp>

#include "coord_converter.h"
#include "lupnt/core/constants.h"
#include "orbit_state.h"

#pragma once

namespace ad = autodiff;

namespace lupnt {

// COE <-> Cart
CartesianOrbitState CoeToCart(const ClassicalOE coe, double mu);
// ad::Vector6real CoeToCart(const ad::Vector6real &coeVec, double mu);
template <typename T>
Eigen::Matrix<T, 6, 1, 0, 6, 1> CoeToCart(
    const Eigen::Matrix<T, 6, 1, 0, 6, 1> &coeVec, double mu);

ClassicalOE CartToCoe(const CartesianOrbitState cartOrbitState, double mu);
ad::Vector6real CartToCoe(const ad::Vector6real &cartVec, double mu);

// ROE <-> COE
ClassicalOE RoeToCoe(const ClassicalOE coeChief, const QuasiNonsingularROE roe);
ad::Vector6real RoeToCoe(const ad::Vector6real &coeChiefVec,
                         const ad::Vector6real &roeVec);

// Inertial <-> RTN
CartesianOrbitState InertialToRtn(
    const CartesianOrbitState &rtnOrigin,
    const CartesianOrbitState &inertialOrbitState);
ad::Vector6real InertialToRtn(const ad::Vector6real &rtnOrigin,
                              const ad::Vector6real &inertialVec);

// COE <-> RTN
CartesianOrbitState CoeToRtn(const ClassicalOE &coeChief,
                             const ClassicalOE &coeDeputy, double mu);
ad::Vector6real CoeToRtn(const ad::Vector6real &coeChiefVec,
                         const ad::Vector6real &coeDeputyVec, double mu);

// COE <-> QNSOE
QuasiNonsingularOE CoeToQnsoe(const ClassicalOE &coe);
ad::Vector6real CoeToQnsoe(const ad::Vector6real &coeVec);

// COE <-> QNSROE
QuasiNonsingularROE QnsoeToQnsroe(const QuasiNonsingularOE &qnsoeChief,
                                  const QuasiNonsingularOE &qnsoeDeputy);
QuasiNonsingularROE CoeToQnsroe(const ClassicalOE &coeChief,
                                const ClassicalOE &coeDeputy);

ClassicalOE QnsoeToCoe(const QuasiNonsingularOE &qnsoe);
ad::Vector6real QnsoeToCoe(const ad::Vector6real &qnsoeVec);

ad::Vector6real EquioeToCoe(const ad::Vector6real &equioe);
ad::Vector6real CoeToEquioe(const ad::Vector6real &coe);

ad::Vector6real MeanToOsculating(const ad::Vector6real &meanCoe, double J2);
ClassicalOE MeanToOsculating(const ClassicalOE &meanCoe, double J2);

ad::Vector6real OsculatingToMean(const ad::Vector6real &osculatingCoe,
                                 double J2);
ClassicalOE OsculatingToMean(const ClassicalOE &osculatingCoe, double J2);

ad::Vector6real CartToQnsoe(const ad::Vector6real &cart, double mu);
QuasiNonsingularOE CartToQnsoe(const CartesianOrbitState &cart, double mu);

QuasiNonsingularOE OscQnsoeToMeanQnsoe(const QuasiNonsingularOE &oscQnsoe,
                                       double J2);

ad::Vector6real DelaunayToCoe(const ad::Vector6real &delaunay, double mu,
                              double n, double t);
ad::Vector6real CoeToDelaunay(const ad::Vector6real &coe, double mu, double n,
                              double t);

std::array<double, 6> ComputeSecondOrderShortPeriod(ad::Vector6real &coe,
                                                    ad::Vector6real &doe);
std::array<double, 6> ComputeFirstOrderMediumPeriod(ad::Vector6real &coe,
                                                    ad::Vector6real &doe);
std::array<double, 6> ComputeSecondOrderMediumPeriod(ad::Vector6real &coe,
                                                     ad::Vector6real &doe);
std::array<double, 6> ComputeCorrectionMediumPeriod(ad::Vector6real &coe,
                                                    ad::Vector6real &doe);

// Anomaly conversions
template <typename T>
T EccentricAnomToTrueAnom(T E, T e);
template <typename T>
T EccentricAnomToMeanAnom(T E, T e);
template <typename T>
T MeanAnomToEccentricAnom(T M, T e);
template <typename T>
T MeanAnomToTrueAnom(T M, T e);
template <typename T>
T TrueAnomToEccentricAnom(T nu, T e);
template <typename T>
T TrueAnomToMeanAnom(T f, T e);

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
    const std::shared_ptr<CartesianOrbitState> state_in, const ad::real epoch,
    const CoordSystem coord_sys_out) {
  auto rv_in = state_in->GetVector();
  auto coord_sys_in = state_in->GetCoordSystem();
  auto rv_out =
      CoordConverter::Convert(rv_in, epoch, coord_sys_in, coord_sys_out);
  auto state_out = std::make_shared<CartesianOrbitState>(rv_out, coord_sys_out);
  return state_out;
}

}  // namespace lupnt