/**
 * @file RadioMeasurement.cpp
 * @author Stanford NAV LAB
 * @brief Class for Radionavigation measurements
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coord_converter.h>

namespace lupnt {

class RadioMeasurement {
 public:
  static real ComputePseudorange(VectorXreal r_tx, VectorXreal r_rx,
                                     real dt_tx, real dt_rx,
                                     real offset);
  static real ComputePseudorangerate(VectorXreal r_tx,
                                         VectorXreal r_rx,
                                         VectorXreal v_tx,
                                         VectorXreal v_rx,
                                         real dt_tx_dot, real dt_rx_dot,
                                         real offset);

  static real ComputeDopplerShift(VectorXreal r_tx,
                                      VectorXreal r_rx,
                                      VectorXreal v_tx,
                                      VectorXreal v_rx, real dt_tx_dot,
                                      real dt_rx_dot, real f,
                                      real offset);
};

}  // namespace lupnt