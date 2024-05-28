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

#include <lupnt/agents/agent.h>
#include <lupnt/core/constants.h>

#include <vector>

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
  VectorXd vis_body;  // visibility body
  VectorXd vis_anttena;

  // Link Budget
  double CN0;

  // delay in two-way ranging [s]
  double delay_tw;

 public:
  ISLMeasurement(Spacecraft* transmit_sat, Spacecraft* receiver_sat);
  double GetCN0() const { return CN0; }

  // True Measurement Generation
  VectorXd GetTrueISLMeasurement(double epoch);

  // Predicted Measurement
  VectorXd GetPredictedISLMeasurement(double epoch, Vector6 rv_trans_pred,
                                      Vector6 rv_rec_pred,
                                      Vector2 clk_trans_pred,
                                      Vector2 clk_rec_pred,
                                      std::vector<ISLMeasurementType> meas_type,
                                      Frame frame_in = Frame::MI);
  VectorXd GetPredictedISLPR(double epoch, Vector6);
};

}  // namespace lupnt