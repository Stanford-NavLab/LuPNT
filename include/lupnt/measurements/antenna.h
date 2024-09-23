/**
 * @file antenna.h
 * @author Stanford NAV LAB
 * @brief Antenna class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <string>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/vector_macros.h"

namespace lupnt {
  class Antenna {
  public:
    Antenna() = default;
    Antenna(const std::string &name) : name_(name) { LoadAntennaPattern(); };

    void LoadAntennaPattern();

    Real ComputeGain(Real theta, Real phi);
    VEC_DEF_REAL_REAL(ComputeGain)

    std::string GetName() { return name_; }
    Real ComputeGain(Vec3 direction);
    MatXd GetGainMatrix() { return gain_; }
    VecXd GetPhiVector() { return phi_; }
    VecXd GetThetaVector() { return theta_; }

  private:
    int n_dim_;         // Number of dimensions (1 or 2)
    std::string name_;  // Name (e.g., Block-IIR_ACE)
    double phi_max_;    // Minimum phiation angle [deg]
    MatXd gain_;        // Gain pattern [dB]
    VecXd phi_;         // Phi angles [deg]
    VecXd theta_;       // Theta angles [deg]

    void FormatAntennaPattern(std::vector<double> &phi, std::vector<double> &theta,
                              std::vector<std::vector<double>> &gain);
  };

}  // namespace lupnt
