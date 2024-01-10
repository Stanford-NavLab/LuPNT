/**
 * @file example_receiver_noise.cc
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-01-08
 *
 * @copyright Copyright (c) 2024
 *
 */

// lupnt includes
#include "example_utils.cc"

using namespace lupnt;
namespace sp = SpiceInterface;

// Util Functions

int main() {
  // setup gps receiver parameters
  GnssReceiverParam rx_param;
  rx_param.Bp = 20.0;
  rx_param.T = 1.0;
  rx_param.Bfe = 20.0;
  rx_param.Bn = 2.0;
  rx_param.D = 0.5;
  rx_param.Rc = 1.023e6;
  rx_param.Tc = 1 / rx_param.Rc;

  // setup gnss measurement
  std::vector<Transmission> trans(0);
  GnssMeasurement gnss_meas(trans);
  gnss_meas.SetGnssReceiverParam(rx_param);

  // Setup C/N0 array
  double CN0_dB = 40.0;

  double pr = gnss_meas.ComputePseudorangeNoise(CN0_dB);
  double prr = gnss_meas.ComputePseudorangeRateNoise(CN0_dB);
  double cp = gnss_meas.ComputeCarrierPhaseNoise(CN0_dB);

  // Print results
  std::cout << "Pseudorange noise at 40 dB-Hz: " << pr << std::endl;
  std::cout << "Pseudorange rate noise at 40 dB-Hz: " << prr << std::endl;
  std::cout << "Carrier phase noise at 40 dB-Hz: " << cp << std::endl;
}