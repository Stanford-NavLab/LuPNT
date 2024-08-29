/**
 * @file gravity_field.h
 * @author Stanford NAVLAB
 * @brief Gravity field model and acceleration calculation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

// C++ includes
#include <string>
#include <vector>

// lupnt includes
#include "lupnt/core/constants.h"
#include "lupnt/core/user_file_path.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"

template <typename T> using Vector3 = Eigen::Vector3<T>;
template <typename T> using Matrix3 = Eigen::Matrix<T, 3, 3>;
template <typename T> using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace lupnt {

  template <typename T> Vector3<T> AccelarationGravityField(const Vector3<T>& r, T GM, T R_ref,
                                                            const MatrixX<T>& CS, int n_max,
                                                            int m_max);

  Vec3 AccelerationPointMass(const Vec3& r, const Vec3& s, Real GM);
  Vec3 AccelerationSolarRadiation(const Vec3& r, const Vec3& r_sun, Real area, Real mass, Real CR,
                                  Real P0, Real AU);
  Vec3 AccelerationDrag(Real mjd_tt, const Vec6& rv, const Mat3& T, Real area, Real mass, Real CD);

  Real Illumination(const Vec3& r, const Vec3& r_sun, Real R_body);
  Real DensityHarrisPriester(Real mjd_tt, const Vec3& r_tod);
  Vec3 AccelerationEarthSpacecraft(Real mjd_tt, const Vec6& rv, Real area, Real mass, Real CR,
                                   Real CD, GravityField<double> grav);
}  // namespace lupnt
