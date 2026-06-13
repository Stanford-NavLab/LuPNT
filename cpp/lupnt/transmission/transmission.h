#pragma once

#include "lupnt/states/state.h"

namespace lupnt {

  struct Transmission {};
  struct Data {};

  struct SimpleGnssData : public Data {
    Cart6 rv_tx;         // [m, m/s] Position and velocity of the transmitter
    ClockState2 clk_tx;  // [s] Clock state of the transmitter
  };

  struct SimpleGnssTransmission : public Transmission {
    Ptr<SimpleGnssData> data;
    Real P_tx;  // [dB-W] Transmitter power
  };

}  // namespace lupnt
