#pragma once

#include "intersatellite_link.h"

#include "radio_measurement.h"

namespace lupnt {

VecXd IslMeasurement::GetTrueISLMeasurement(double epoch) {
  VecXd isl_meas;
  return isl_meas;
}

VecXd IslMeasurement::GetPredictedISLMeasurement(
    double epoch, Vec6 rv_trans_pred, Vec6 rv_rec_pred, Vec2 clk_trans_pred,
    Vec2 clk_rec_pred, std::vector<ISLMeasurementType> meas_type,
    Frame frame_in = Frame::MOON_CI) {
  // Measurement vector
  VecX z(n_meas * meas_type.size());
  int idx = 0;

  for (auto type : meas_type) {
    switch (type) {
      case IslMeasurementType::PR:
        z.segment(idx, n_meas) =
            ComputePseudorange(rv_trans_pred.head(3), rv_rec_pred.head(3),
                               clk_trans_pred(0), clk_rec_pred(0), 0.0);
        idx += 1;
        break;
      case IslMeasurementType::PRR:
        z.segment(idx, n_meas) = GetPseudorangerate(with_noise, seed);
        break;
      case IslMeasurementType::TWR:
        z.segment(idx, n_meas) = GetCarrierPhase(with_noise, seed);
        break;
      case IslMeasurementType::TWRR:
        z.segment(idx, n_meas) = GetPseudorange(with_noise, seed);
        break;
      case IslMeasurementType::DWR:
        z.segment(idx, n_meas) = GetPseudorangerate(with_noise, seed);
        break;
      case IslMeasurementType::DWRR:
        z.segment(idx, n_meas) = GetCarrierPhase(with_noise, seed);
        break;
      default:
        std::cout << "Measurement type: " << type << " not supported"
                  << std::endl;
        break;
    }
  }
  return z;

}  // namespace lupnt

}  // namespace lupnt