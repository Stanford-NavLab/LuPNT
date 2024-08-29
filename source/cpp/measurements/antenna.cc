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

  /**
   * @brief Set the Antenna Pattern object
   *
   * GPS SV(transmitter) antenna file definitions:
   * The format for 1D (elevation) antenna patttern files is as follows:
   * % angle [deg] gain [db]
   * 0.0  gain1
   * . .
   * . .
   * angN gainN
   *
   * The format for 2D antenna patttern files is as follows: (angles in deg, gain
   * in dB)
   *
   * -50  0.0  az2  ...  azM
   * 0.0  gain1,1 gain1,2 ...  gain1,M
   * . .
   * . .
   * elN  gainN,1 gainN,2 ...  gainN,M
   *
   */
  void Antenna::LoadAntennaPattern() {
    // Check for omni-directional antenna
    if (comms_name_ == "") return;

    // Find file
    std::filesystem::path antennaPath(GetDataPath() / "antenna");
    auto filePath = FindFileInDir(antennaPath, comms_name_ + ".txt");
    assert(filePath.has_value() && "Antenna pattern file not found.");

    // Open file
    std::ifstream file(filePath->string(), std::ifstream::in);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open antenna pattern file: " + filePath->string());
    }

    // Read the header
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] != '%') {
        break;
      }
    }

    // Read the pattern
    /**
     * 2D antenna pattern
     * -50  0.0     az2     ...  azM
     * 0.0  gain1,1 gain1,2 ...  gain1,M
     * . .
     * . .
     * elN  gainN,1 gainN,2 ...  gainN,M
     *
     * 1D antenna pattern
     * 0.0  gain1
     * . .
     * . .
     * angN gainN
     * */

    // Read the first line
    std::stringstream ss(line);
    std::string token;
    std::vector<double> tokens;
    while (std::getline(ss, token, ',')) {
      tokens.push_back(std::stod(token));
    }

    // Check if 1D or 2D
    if (tokens.size() == 2) {
      // 1D
      // Read the rest of the file
      std::vector<double> angles;
      std::vector<double> gains;
      while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<double> tokens;
        while (std::getline(ss, token, ',')) {
          tokens.push_back(std::stod(token));
        }
        angles.push_back(tokens[0]);
        gains.push_back(tokens[1]);
      }

      // Create the pattern
      antenna_pattern_.resize(2, angles.size());
      for (size_t i = 0; i < angles.size(); i++) {
        antenna_pattern_(0, i) = angles[i];
        antenna_pattern_(1, i) = gains[i];
      }

    } else {
      // 2D
      // Read the rest of the file
      std::vector<double> thetas;
      std::vector<double> phis;
      std::vector<std::vector<double>> gains;

      for (size_t i = 1; i < tokens.size(); i++) {
        phis.push_back(tokens[i]);
      }
      // Add the first phi again for 0-360 deg
      phis.push_back(360);
      while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<double> tokens;
        while (std::getline(ss, token, ',')) {
          tokens.push_back(std::stod(token));
        }
        thetas.push_back(tokens[0]);
        gains.push_back(std::vector<double>(tokens.begin() + 1, tokens.end()));
      }

      // Create the pattern
      antenna_pattern_.resize(int(thetas.size()) + 1, int(phis.size()) + 1);
      antenna_pattern_(0, 0) = -50;
      for (size_t i = 0; i < phis.size(); i++) {
        antenna_pattern_(0, i + 1) = phis[i];
      }
      // Add the first phi again for 0-360 deg
      // antenna_pattern_(0, phis.size() + 1) = phis[0];

      for (size_t i = 0; i < thetas.size(); i++) {
        antenna_pattern_(i + 1, 0) = thetas[i];
        for (size_t j = 0; j < phis.size() - 1; j++) {
          antenna_pattern_(i + 1, j + 1) = gains[i][j];
        }
        // Add the first phi again for 360 deg
        antenna_pattern_(i + 1, phis.size()) = gains[i][0];
      }
      file.close();
    }
  }

  /**
   * @brief Get the Antenna Gain object
   *
   * @param theta elevation angle [deg]
   * @param phi azimuth angle [deg]
   * @return double   antenna gain [dB]
   */
  double Antenna::GetAntennaGain(double theta, double phi) {
    double gain = 0.0;

    // Check omni-directional antennas
    if (antenna_pattern_.size() == 0) return gain;

    // The anglemask is the max theta angle for which gain will be evaluated. It
    // is the minimum of the max defined angle for the antenna antenna_pattern_.
    double anglemask = antenna_pattern_.col(0).maxCoeff();
    if (theta > anglemask) return -500;

    // Determine if the pattern is elevation only (1-D) or azimuth and elevation
    // (2-D)
    if (antenna_pattern_.cols() > 2) {
      // Use both transmitter azimuth and elevation to compute gain from a 2D
      // antenna model
      VecXd x = antenna_pattern_.col(0).segment(1, antenna_pattern_.rows() - 1);
      VecXd y = antenna_pattern_.row(0).segment(1, antenna_pattern_.cols() - 1);
      MatXd z
          = antenna_pattern_.block(1, 1, antenna_pattern_.rows() - 1, antenna_pattern_.cols() - 1);
      double res = LinearInterp2d(x, y, z, theta, phi + 180.0);
      gain = res;
    } else if (antenna_pattern_.cols() == 2) {
      // Use only transmitter elevation angle to compute gain from a 1D antenna
      // model
      VecXd x = antenna_pattern_.col(0);
      VecXd y = antenna_pattern_.col(1);
      double res = LinearInterp1d(x, y, theta);
      gain = res;
    }

    return gain;
  }

  double Antenna::GetAntennaGain(Vec3d direction) {
    direction.normalize();

    // Compute elevation and azimuth angles
    double theta = acos(direction.dot(Vec3d::UnitZ())) * DEG;
    double phi = atan2(direction.y(), direction.x()) * DEG;

    return GetAntennaGain(theta, phi);
  }

}  // namespace lupnt
