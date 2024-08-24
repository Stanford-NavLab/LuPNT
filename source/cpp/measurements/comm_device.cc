/**
 * @file comm_device.cc
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-08-20
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/measurements/comm_device.h"

#include "lupnt/measurements/radio_measurement.h"

namespace lupnt {
std::shared_ptr<Transmitter> CreateTransmitter(const std::string name, double P_tx, double freq,
                              double bandwidth) {

  std::shared_ptr<Transmitter> tx;
  tx->name = name;
  tx->freq_tx = freq;
  tx->P_tx = P_tx;
  tx->bandwidth = bandwidth;

  return tx;
}



}  // namespace lupnt
