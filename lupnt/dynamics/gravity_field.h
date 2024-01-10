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

std::tuple<double, double> spharm_vwmm(int m_in, double Vm_1m_1, double Wm_1m_1,
                                       const Vector3d &x_R, double Re);

std::tuple<double, double> spharm_vwm1m(int m_in, double Vmm, double Wmm,
                                        const Vector3d &x_R, double Re);

std::tuple<double, double> spharm_vwnm(int n_in, int m_in, double Vn_1m,
                                       double Vn_2m, double Wn_1m, double Wn_2m,
                                       const Vector3d &x_R, double Re);

Vector3d Facc_j(const Vector3d &facc_R, const Matrix3d &Ur2j);

std::tuple<Vector3d, Vector3d> spharm_dvwdx(int n_in, int m_in, double Vn1m,
                                            double Vn1m1, double Vn1m_1,
                                            double Wn1m, double Wn1m1,
                                            double Wn1m_1, double Re);

std::tuple<Matrix3d, Matrix3d> spharm_d2vwdx2(int n_in, int m_in, double Vn2m,
                                              double Vn2m1, double Vn2m2,
                                              double Vn2m_1, double Vn2m_2,
                                              double Wn2m, double Wn2m1,
                                              double Wn2m2, double Wn2m_1,
                                              double Wn2m_2, double Re);

Vector3d spharm_acc_ecr(int nmax, int mmax, const Vector3 &x_R_in, double Re,
                        double GMe, const MatrixXd &Cnm, const MatrixXd &Snm);

}  // namespace lupnt