#include "lupnt/measurements/isl_measurement.h"

#include "lupnt/measurements/radio_measurement.h"

namespace lupnt {

VecX IslMeasurement::GetTrueIslMeasurement(double epoch) {
  VecX isl_meas;
  return isl_meas;
}

VecX IslMeasurement::GetPredictedIslMeasurement(
    double epoch, Vec6 rv_trans_pred, Vec6 rv_rec_pred, Vec2 clk_trans_pred,
    Vec2 clk_rec_pred, std::vector<IslMeasurementType> meas_type,
    Frame frame_in) {
  // Measurement vector
  VecX z(meas_type.size());
  int idx = 0;
  bool with_noise = false;
  int seed = 0;

  // for (auto type : meas_type) {
  //   switch (type) {
  //     case IslMeasurementType::PR:
  //       z.segment(idx, n_meas) =
  //           ComputePseudorange(rv_trans_pred.head(3), rv_rec_pred.head(3),
  //                              clk_trans_pred(0), clk_rec_pred(0), 0.0);
  //       idx += 1;
  //       break;
  //     case IslMeasurementType::PRR:
  //       z.segment(idx, n_meas) = GetPseudorangerate(with_noise, seed);
  //       break;
  //     case IslMeasurementType::TWR:
  //       z.segment(idx, n_meas) = GetCarrierPhase(with_noise, seed);
  //       break;
  //     case IslMeasurementType::TWRR:
  //       z.segment(idx, n_meas) = GetPseudorange(with_noise, seed);
  //       break;
  //     case IslMeasurementType::DWR:
  //       z.segment(idx, n_meas) = GetPseudorangerate(with_noise, seed);
  //       break;
  //     case IslMeasurementType::DWRR:
  //       z.segment(idx, n_meas) = GetCarrierPhase(with_noise, seed);
  //       break;
  //     default:
  //       std::cout << "Measurement type: " << type << " not supported"
  //                 << std::endl;
  //       break;
  //   }
  // }

  return z;
}

}  // namespace lupnt