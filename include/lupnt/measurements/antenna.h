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

namespace lupnt {
  class Antenna {
  public:
    std::string comms_name_;            // Name of the antenna
    MatXd antenna_pattern_;             // Antenna gain pattern [deg & dB]
    double antenna_mask_ = 80.0 * RAD;  // Cut off angle for the transmit antenna [rad]

    Antenna() = default;
    Antenna(std::string comms_name) : comms_name_(comms_name) { LoadAntennaPattern(); };

    void LoadAntennaPattern();
    double GetAntennaGain(double theta, double phi);
    double GetAntennaGain(Vec3d direction);
  };

}  // namespace lupnt
