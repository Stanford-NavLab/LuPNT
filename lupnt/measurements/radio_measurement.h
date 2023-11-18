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
  static ad::real ComputePseudorange(ad::VectorXreal r_tx, ad::VectorXreal r_rx,
                                     ad::real dt_tx, ad::real dt_rx,
                                     ad::real offset);
  static ad::real ComputePseudorangerate(ad::VectorXreal r_tx,
                                         ad::VectorXreal r_rx,
                                         ad::VectorXreal v_tx,
                                         ad::VectorXreal v_rx,
                                         ad::real dt_tx_dot, ad::real dt_rx_dot,
                                         ad::real offset);

  static ad::real ComputeDopplerShift(ad::VectorXreal r_tx,
                                      ad::VectorXreal r_rx,
                                      ad::VectorXreal v_tx,
                                      ad::VectorXreal v_rx, ad::real dt_tx_dot,
                                      ad::real dt_rx_dot, ad::real f,
                                      ad::real offset);
};

}  // namespace lupnt