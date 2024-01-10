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
  rx_param.Bp = 5.0;
  rx_param.T = 20e-3;
  rx_param.b = 2;  // normalized bandwidth
  rx_param.Bn = 0.2;
  rx_param.D = 0.3;

  // lambda
  double freq = 1575.42e6;
  double lambda = C / freq;

  // setup gnss measurement
  std::vector<Transmission> trans(1);
  trans[0].chip_rate = 1.023e6;  // L1 chip rate

  GnssMeasurement gnss_meas(trans);
  gnss_meas.SetGnssReceiverParam(rx_param);

  // Setup C/N0 array
  std::vector<double> CN0_dB_vec{15.0, 20.0, 30.0, 40.0};

  // Print results
  for (auto CN0_dB : CN0_dB_vec) {
    double pr = gnss_meas.ComputePseudorangeNoise(CN0_dB);
    double prr = gnss_meas.ComputePseudorangeRateNoise(CN0_dB, lambda);
    double cp = gnss_meas.ComputeCarrierPhaseNoise(CN0_dB, lambda);

    // Print results
    std::cout << "C/N0 [dB-Hz] : " << CN0_dB << std::endl;
    std::cout << "  Pseudorange noise      [m]   : " << 1000 * pr << std::endl;
    std::cout << "  Pseudorange rate noise [m/s] : " << 1000 * prr << std::endl;
    std::cout << "  Carrier phase noise    [m]   : " << 1000 * cp << std::endl;
    std::cout << " " << std::endl;
  }
}