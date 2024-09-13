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
      std::vector<double> elev, gain;
      while (std::getline(file, line)) {
        ss = std::stringstream(line);
        tokens.clear();
        while (std::getline(ss, token, ',')) tokens.push_back(std::stod(token));
        elev.push_back(tokens[0]);
        gain.push_back(tokens[1]);
      }
      gain_ = MatXd::Zero(elev.size(), 1);
      elev_ = VecXd::Zero(elev.size());
      for (size_t i = 0; i < gain.size(); i++) {
        gain_(i, 0) = gain[i];
        elev_(i) = elev[i];
      }
    } else {
      // 2D
      n_dim_ = 2;
      std::vector<double> elev, azim;
      std::vector<std::vector<double>> gain;
      for (size_t i = 1; i < tokens.size(); i++) azim.push_back(tokens[i]);
      while (std::getline(file, line)) {
        ss = std::stringstream(line);
        tokens.clear();
        while (std::getline(ss, token, ',')) tokens.push_back(std::stod(token));
        elev.push_back(tokens[0]);
        gain.push_back(std::vector<double>(tokens.begin() + 1, tokens.end()));
      }
      azim.push_back(360.0);
      for (size_t i = 0; i < gain.size(); i++) gain[i].push_back(gain[i][0]);

      gain_ = MatXd::Zero(elev.size(), azim.size());
      elev_ = VecXd::Zero(elev.size());
      azim_ = VecXd::Zero(azim.size());
      for (size_t i = 0; i < gain.size(); i++) {
        elev_(i) = elev[i];
        for (size_t j = 0; j < gain[i].size(); j++) gain_(i, j) = gain[i][j];
      }
      for (size_t i = 0; i < azim.size(); i++) azim_(i) = azim[i];
    }
    file.close();
  }

  /**
   * @brief Get the Antenna Gain object
   *
   * @param elev elevation angle [rad]
   * @param az azimuth angle [rad]
   * @return double   antenna gain [dB]
   */
  double Antenna::ComputeGain(Real elev, Real azim) {
    elev *= DEG;
    azim *= DEG;
    double gain = 0.0;

    // Omni-directional antenna
    if (n_dim_ == 0) return gain;

    // Maximum elevation angle
    if (elev > elev_.maxCoeff()) return -500;

    // Gain
    if (n_dim_ == 2)
      gain = LinearInterp2d(elev_, azim_, gain_, elev.val(), azim.val());
    else if (n_dim_ == 1)
      gain = LinearInterp1d(elev_, gain_, elev.val());

    return gain;
  }

  VEC_IMP_REAL_REAL(ComputeGain)

  double Antenna::ComputeGain(Vec3 direction) {
    direction.normalize();

    // Compute elevation and azimuth angles
    Real elev = acos(direction.dot(Vec3::UnitZ()));
    Real azim = atan2(direction.y(), direction.x());

    return ComputeGain(elev, azim);
  }

}  // namespace lupnt
