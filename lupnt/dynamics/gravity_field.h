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

std::tuple<Real, Real> spharm_vwmm(int m_in, Real Vm_1m_1, Real Wm_1m_1,
                                   const Vec3 &x_R, double Re);

std::tuple<Real, Real> spharm_vwm1m(int m_in, Real Vmm, Real Wmm,
                                    const Vec3 &x_R, double Re);

std::tuple<Real, Real> spharm_vwnm(int n_in, int m_in, Real Vn_1m, Real Vn_2m,
                                   Real Wn_1m, Real Wn_2m, const Vec3 &x_R,
                                   double Re);

Vec3 Facc_j(const Vec3 &facc_R, const Mat3 &Ur2j);

std::tuple<Vec3, Vec3> spharm_dvwdx(int n_in, int m_in, Real Vn1m, Real Vn1m1,
                                    Real Vn1m_1, Real Wn1m, Real Wn1m1,
                                    Real Wn1m_1, double Re);

std::tuple<Mat3d, Mat3d> spharm_d2vwdx2(int n_in, int m_in, double Vn2m,
                                        double Vn2m1, double Vn2m2,
                                        double Vn2m_1, double Vn2m_2,
                                        double Wn2m, double Wn2m1, double Wn2m2,
                                        double Wn2m_1, double Wn2m_2,
                                        double Re);

Vec3 spharm_acc_ecr(int nmax, int mmax, const Vec3 &x_R_in, double Re,
                    double GMe, const VecXd &Cnm, const VecXd &Snm);

}  // namespace lupnt