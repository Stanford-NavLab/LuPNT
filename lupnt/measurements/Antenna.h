
#pragma once

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <string>

#include "lupnt/core/Constants.h"

namespace LPT {
class Antenna {
 public:
  std::string comms_name_;           // Name of the antenna
  Eigen::MatrixXd antenna_pattern_;  // Antenna gain pattern [deg & dB]
  double antenna_mask_ =
      80.0 * RAD_PER_DEG;  // Cut off angle for the transmit antenna [rad]

  Antenna() = default;
  Antenna(std::string comms_name) : comms_name_(comms_name) {
    LoadAntennaPattern();
  };

  void LoadAntennaPattern();
  double GetAntennaGain(double theta, double phi);
  double GetAntennaGain(Eigen::Vector3d direction);
};

}  // namespace LPT