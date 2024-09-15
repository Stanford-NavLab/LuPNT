/**
 * @file Antenna.cpp
 * @author Stanford NAV LAB
 * @brief Antenna class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "lupnt/measurements/antenna.h"

#include <filesystem>
#include <iostream>
#include <sstream>

#include "lupnt/core/file.h"
#include "lupnt/numerics/interpolation.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  void Antenna::FormatAntennaPattern(std::vector<double> &phi, std::vector<double> &theta,
                                     std::vector<std::vector<double>> &gain) {
    // Format theta to [0, 360] and phi to [-180, 180]
    if (theta.front() == 0 && theta.back() <= 360 && phi.back() <= 180) {
      if (phi.back() >= 80 && phi.back() < 90) {
        phi.push_back(90);
        gain.push_back(gain.front());
      }
      if (phi.back() >= 170 && phi.back() < 180) {
        phi.push_back(180);
        gain.push_back(gain.front());
      }
      phi_max_ = phi.back();
      if (theta.back() < 360) {
        theta.push_back(360);
        for (auto &g : gain) g.push_back(g.front());
      }
    } else if (theta.front() == -180 && theta.back() == 180 && phi.front() == 0
               && phi.back() == 180) {
      std::vector<double> theta_tmp, g_tmp;
      std::vector<std::vector<double>> gain_tmp;
      for (size_t i = 0; i < phi.size(); i++) {
        g_tmp.clear();
        for (size_t j = 0; j < theta.size(); j++) {
          if (theta[j] >= 0) {  // [0, 180]
            if (i == 0) theta_tmp.push_back(theta[j]);
            g_tmp.push_back(gain[i][j]);
          }
        }
        for (size_t j = 1; j < theta.size(); j++) {
          if (theta[j] <= 0) {  // (180, 360]
            if (i == 0) theta_tmp.push_back(theta[j] + 360);
            g_tmp.push_back(gain[i][j]);
          }
        }
        gain_tmp.push_back(g_tmp);
      }
      phi_max_ = phi.back();
      theta = theta_tmp;
      gain = gain_tmp;
    } else {
      throw std::runtime_error("Invalid antenna pattern");
    }
  }
  void Antenna::LoadAntennaPattern() {
    // Check for omni-directional antenna
    if (name_ == "") {
      n_dim_ = 0;
      return;
    }

    // File
    std::filesystem::path filePath = GetFilePath(name_);
    std::ifstream file(filePath, std::ifstream::in);
    assert(file.is_open() && "File not found.");

    // Header
    std::string line;
    while (std::getline(file, line))
      if (line[0] != '%') break;

    // 2D antenna pattern
    // ----------------------------------
    // -50  0.0     az2     ...  azM
    // 0.0  gain1,1 gain1,2 ...  gain1,M
    // . .
    // . .
    // elN  gainN,1 gainN,2 ...  gainN,M
    //
    // 1D antenna pattern
    // ----------------------------------
    // 0.0  gain1
    // . .
    // . .
    // angN gainN

    // Read the first line
    std::stringstream ss(line);
    std::string token;
    std::vector<double> tokens;
    while (std::getline(ss, token, ',')) tokens.push_back(std::stod(token));

    // Check if 1D or 2D
    if (tokens.size() == 2) {
      // 1D pattern
      n_dim_ = 1;
      std::vector<double> phi, gain;
      while (std::getline(file, line)) {
        ss = std::stringstream(line);
        tokens.clear();
        while (std::getline(ss, token, ',')) tokens.push_back(std::stod(token));
        phi.push_back(tokens[0]);
        gain.push_back(tokens[1]);
      }
      gain_ = MatXd::Zero(phi.size(), 1);
      phi_ = VecXd::Zero(phi.size());
      for (size_t i = 0; i < gain.size(); i++) {
        gain_(i, 0) = gain[i];
        phi_(i) = Wrap2Pi(phi[i] * RAD).val();
      }
    } else {
      // 2D
      n_dim_ = 2;
      std::vector<double> phi, theta;
      std::vector<std::vector<double>> gain;
      for (size_t i = 1; i < tokens.size(); i++) theta.push_back(tokens[i]);
      while (std::getline(file, line)) {
        ss = std::stringstream(line);
        tokens.clear();
        while (std::getline(ss, token, ',')) tokens.push_back(std::stod(token));
        phi.push_back(tokens[0]);
        gain.push_back(std::vector<double>(tokens.begin() + 1, tokens.end()));
      }

      // Format theta to [0, 360] and phi to [-180, 180]
      FormatAntennaPattern(phi, theta, gain);

      // Pattern
      gain_ = MatXd::Zero(phi.size(), theta.size());
      phi_ = VecXd::Zero(phi.size());
      theta_ = VecXd::Zero(theta.size());
      for (size_t i = 0; i < gain.size(); i++) {
        phi_(i) = phi[i];
        for (size_t j = 0; j < gain[i].size(); j++) gain_(i, j) = gain[i][j];
      }
      for (size_t i = 0; i < theta.size(); i++) theta_(i) = theta[i];
    }
    file.close();
  }

  Real Antenna::ComputeGain(Real theta, Real phi) {
    Real gain = 0.0;
    if (n_dim_ == 0) return gain;  // Omni-directional antenna

    theta = Wrap2TwoPi(theta) * DEG;
    phi = Wrap2Pi(phi) * DEG;
    if (phi > phi_max_ || phi < -phi_max_) return NAN;  // Minimum angle
    if (phi < 0 && phi_(0) >= 0) phi = -phi;            // Symmetry

    // Interpolation
    if (n_dim_ == 2) {
      gain = LinearInterp2d(phi_, theta_, gain_, phi.val(), theta.val());
    } else if (n_dim_ == 1) {
      gain = LinearInterp1d(phi_, gain_, phi.val());
    }
    return gain;
  }

  VEC_IMP_REAL_REAL(Antenna::ComputeGain)

  Real Antenna::ComputeGain(Vec3 direction) {
    direction.normalize();

    // Compute phiation and thetauth angles
    Real phi = acos(direction.dot(Vec3::UnitZ()));
    Real theta = atan2(direction.y(), direction.x());

    return ComputeGain(theta, phi);
  }

}  // namespace lupnt
