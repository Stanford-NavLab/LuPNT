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

#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  class RadioMeasurement {
  public:
    static Real ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx, Real dt_rx, Real offset);

    static Real ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                                       Real dt_rx_dot, Real offset);

    static Real ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                                    Real dt_rx_dot, Real f, Real offset);
  };

}  // namespace lupnt
