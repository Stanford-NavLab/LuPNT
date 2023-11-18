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
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// autodiff includes
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

// lupnt includes
#include "lupnt/core/user_file_path.h"
#include "lupnt/numerics/math_utils.h"

namespace ad = autodiff;

namespace lupnt {

enum class BodyId {
  VENUS = 299,
  MOON = 301,
  EARTH = 399,
  MARS = 499,
};

struct BodyData {
  std::string name;
  double GM;
  double R;

  std::string filepath;
  int headerlines;
  std::string delimiter;
};

std::vector<std::string> split(const std::string &s, char delimiter);

void ReadData(const std::string &filepath, int N, int headerlines,
              const std::string &delimiter, std::vector<int> &idN,
              std::vector<int> &idM, std::vector<double> &C,
              std::vector<double> &S, std::vector<double> &sigC,
              std::vector<double> &sigS);

BodyData GetBodyData(const BodyId bodyID);

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> LoadGravityCoefficients(
    BodyData bd, int nmax);

std::tuple<double, double> spharm_vwmm(int m_in, double Vm_1m_1, double Wm_1m_1,
                                       const Eigen::Vector3d &x_R, double Re);

std::tuple<double, double> spharm_vwm1m(int m_in, double Vmm, double Wmm,
                                        const Eigen::Vector3d &x_R, double Re);

std::tuple<double, double> spharm_vwnm(int n_in, int m_in, double Vn_1m,
                                       double Vn_2m, double Wn_1m, double Wn_2m,
                                       const Eigen::Vector3d &x_R, double Re);

Eigen::Vector3d Facc_j(const Eigen::Vector3d &facc_R,
                       const Eigen::Matrix3d &Ur2j);

std::tuple<Eigen::Vector3d, Eigen::Vector3d> spharm_dvwdx(
    int n_in, int m_in, double Vn1m, double Vn1m1, double Vn1m_1, double Wn1m,
    double Wn1m1, double Wn1m_1, double Re);

std::tuple<Eigen::Matrix3d, Eigen::Matrix3d> spharm_d2vwdx2(
    int n_in, int m_in, double Vn2m, double Vn2m1, double Vn2m2, double Vn2m_1,
    double Vn2m_2, double Wn2m, double Wn2m1, double Wn2m2, double Wn2m_1,
    double Wn2m_2, double Re);

Eigen::Vector3d spharm_acc_ecr(int nmax, int mmax,
                               const ad::Vector3real &x_R_in, double Re,
                               double GMe, const Eigen::MatrixXd &Cnm,
                               const Eigen::MatrixXd &Snm);

}  // namespace lupnt