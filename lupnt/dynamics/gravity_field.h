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

std::tuple<real, real> spharm_vwmm(int m_in, real Vm_1m_1, real Wm_1m_1,
                                   const Vector3 &x_R, double Re);

std::tuple<real, real> spharm_vwm1m(int m_in, real Vmm, real Wmm,
                                    const Vector3 &x_R, double Re);

std::tuple<real, real> spharm_vwnm(int n_in, int m_in, real Vn_1m, real Vn_2m,
                                   real Wn_1m, real Wn_2m, const Vector3 &x_R,
                                   double Re);

Vector3 Facc_j(const Vector3 &facc_R, const Matrix3 &Ur2j);

std::tuple<Vector3, Vector3> spharm_dvwdx(int n_in, int m_in, real Vn1m,
                                          real Vn1m1, real Vn1m_1, real Wn1m,
                                          real Wn1m1, real Wn1m_1, double Re);

std::tuple<Matrix3d, Matrix3d> spharm_d2vwdx2(int n_in, int m_in, double Vn2m,
                                              double Vn2m1, double Vn2m2,
                                              double Vn2m_1, double Vn2m_2,
                                              double Wn2m, double Wn2m1,
                                              double Wn2m2, double Wn2m_1,
                                              double Wn2m_2, double Re);

Vector3 spharm_acc_ecr(int nmax, int mmax, const Vector3 &x_R_in, double Re,
                       double GMe, const MatrixXd &Cnm, const MatrixXd &Snm);

}  // namespace lupnt