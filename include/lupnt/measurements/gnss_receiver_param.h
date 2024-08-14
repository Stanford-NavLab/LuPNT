/**
 * @file gnss_receiver_param.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-01-05
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

namespace lupnt {

  struct GnssReceiverParam {
    // Receiver chip parameter
    double Bp;  // Carrier loop noise bandwidth [Hz]
    double T;   // Tracking loop integration time [s]
    double b;   // normalized bandwidth [Hz]
    double Bn;  // Code loop noise bandwidth [Hz]
    double D;   // Early-to-late correlator spacing (chips)
    double Rc;  // chipping rate
    double Tc;  // chip period [s]
  };

}  // namespace lupnt
