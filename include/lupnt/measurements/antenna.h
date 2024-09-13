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
    Antenna(const std::string& name) : name_(name) { LoadAntennaPattern(); };

    void LoadAntennaPattern();

    Real ComputeGain(Real azim, Real elev);
    VEC_DEF_REAL_REAL(ComputeGain)

    Real ComputeGain(Vec3 direction);
    MatXd GetGainPattern() { return gain_; }
    VecXd GetElevationAngles() { return elev_; }
    VecXd GetAzimuthAngles() { return azim_; }

  private:
    int n_dim_;         // Number of dimensions (1 or 2)
    std::string name_;  // Name (e.g., Block-IIR_ACE)
    MatXd gain_;        // Gain pattern [dB]
    VecXd elev_;        // Elevation angles [deg]
    VecXd azim_;        // Azimuth angles [deg]
  };

}  // namespace lupnt
