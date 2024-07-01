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

namespace lupnt {
Vec3 AccelarationGravityField(const Vec3& r, const Mat3& E, Real GM, Real R_ref,
                              const MatX& CS, int n_max, int m_max);
Vec3d AccelarationGravityFieldEigen(const Vec3d& r, const Mat3d& E, double GM,
                                    double R_ref, const MatXd& CS, int n_max,
                                    int m_max);
Vec3 AccelerationPointMass(const Vec3& r, const Vec3& s, double GM);
Vec3 AccelerationSolarRadiation(const Vec3& r, const Vec3& r_sun, Real area,
                                Real mass, Real CR, Real P0, Real AU);
Vec3 AccelerationDrag(Real mjd_tt, const Vec6& rv, const Mat3& T, Real area,
                      Real mass, Real CD);

Real DensityHarrisPriester(Real mjd_tt, const Vec3& r_tod);
Vec3 AccelerationEarthSpacecraft(Real mjd_tt, const Vec6& rv, Real area,
                                 Real mass, Real CR, Real CD,
                                 GravityField grav);
}  // namespace lupnt