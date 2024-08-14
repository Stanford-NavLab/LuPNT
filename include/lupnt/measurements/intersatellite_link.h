/**
 * @file intersatellite_link.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-04-04
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "radio_measurement.h"

namespace lupnt {

  enum class ISLMeasurementType {
    PR,   // Pseudorange
    PRR,  // Pseudorange rate
    TWR,  // Two-way range
    TWRR  // Two-way range rate
  };

  class ISLMeasurement {
  private:
    Spacecraft* transmit_sat;
    Spacecraft* receiver_sat;

    // visibility
    std::vector<NaifId> occult_planets;
    VecXd vis_body;  // visibility body
    VecXd vis_anttena;

    // Link Budget
    double CN0;

    // delay in two-way ranging [s]
    double delay_tw;

  public:
    ISLMeasurement(Spacecraft* transmit_sat, Spacecraft* receiver_sat);
    double GetCN0() const { return CN0; }

    // True Measurement Generation
    VecXd GetTrueISLMeasurement(double epoch);

    // Predicted Measurement
    VecXd GetPredictedISLMeasurement(double epoch, Vec6 rv_trans_pred, Vec6 rv_rec_pred,
                                     Vec2 clk_trans_pred, Vec2 clk_rec_pred,
                                     std::vector<ISLMeasurementType> meas_type,
                                     Frame frame_in = Frame::MOON_CI);
    VecXd GetPredictedISLPR(double epoch, Vec6);
  };

}  // namespace lupnt
