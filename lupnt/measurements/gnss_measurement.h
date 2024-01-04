/**
 * @file gnss_measurement.h
 * @author Stanford NAV LAB
 * @brief Class that constructs GPS measurement data from channel observables
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <lupnt/core/constants.h>
#include <lupnt/measurements/transmission.h>

#include <memory>

namespace lupnt {

class GnssMeasurement {
  // Implemenation based on Gnss SDR Observables block:
  // https://gnss-sdr.org/docs/sp-blocks/observables/
 private:
  int n_meas = 0;                         // Number of measurements
  std::vector<Transmission> trans_store;  // list of transmittion data
  std::vector<int> ID_tx;                 // ID of the transmitter (n_meas)

  // Visbility
  VectorXd vis_earth;
  VectorXd vis_moon;
  VectorXd vis_antenna;
  VectorXd vis_atmos;
  VectorXd vis_ionos;

  double c = 299792.458;  // Speed of light [km/s]

  // Pseudorange measurement
  double t_rx;  // Signal reception time measured by the receiver clock [s]
  double
      t_tx;  // Signal transmission time measured by the transmitter clock [s]

  VectorX P_rx;     // Pseudorange measurement [km] (n_meas * n_bands)
  VectorXd rho_rx;  // True range from the transmitter to the receiver’s
                    // antenna [km] (n_meas)

  double dt_rx;    // Receiver clock offset from Gnss time [s]
  VectorXd dt_tx;  // Transmitter clock offset from Gnss time [s] (n_meas)

  VectorXd I_rx;   // Ionospheric delay [km] (n_meas * n_bands)
  VectorXd T_rx;   // Tropospheric delay [km] (n_meas)
  VectorXd eps_P;  // Pseudorange measurement noise [km] (n_meas)

  // Carrier phase measurement
  VectorX phi_rx;      // Carrier phase measurement [cycles] (n_meas * n_bands)
  VectorXd phi_rx_tx;  // Phase of the receiver's local oscillator at
                       // time t [cycles] (n_meas * n_bands)
  VectorXd phi_tx;     // Phase of the transmitted signal at time t [cycles]
                       // (n_meas * n_bands)
  double t;
  VectorXd N_rx;     // Carrier phase int ambiguity [cycles] (n_meas * n_bands)
  VectorXd eps_phi;  // Carrier phase measurement noise [cycles]
                     // (n_meas * n_bands)
  VectorXd f;  // Carrier frequency [Hz] (n_bands) // TODO: Convert to vector
  double t_0;  // Initial time [s]
  VectorXd phi_0;  // Initial phase of the transmitted signal at time t_0
                   // [cycles] (n_meas * n_bands)
  double
      lambda;  // Carrier wavelength [km] (n_bands) // TODO: Convert to vector
  VectorXd L_rx_txr;  // Integer component of the receiver’s numerically
                      // controlled oscillator (NCO) phase
  VectorXd K_rx_txr;  // Integer component of the propagation term

  // Phase-range masurement
  VectorXd Phi_rx;     // Phase-range measurement [km] (n_meas * n_bands)
                       // (n_meas * n_bands)
  VectorXd B_rx;       // Carrier phase bias [cycles] (n_meas * n_bands)
  VectorXd dPhi_rx;    // Receiver’s antenna phase center variation [km]
                       // (n_meas * n_bands)
  VectorXd d_rx_pco;   // Receiver’s antenna phase center offset in local
                       // coordinates [km] (n_meas) (n_bands x 3)
  VectorXd d_rx_pcv;   // Receiver’s antenna phase center variation [km]
                       // (n_meas * n_bands)
  VectorXd d_tx_pco;   // Transmitter’s antenna phase center offset in
                       // local coordinates [km] (n_meas * n_bands x 3)
  VectorXd d_tx_pcv;   // Transmitter’s antenna phase center variation
                       // [km] (n_meas * n_bands)
  MatrixXd e_rx_enu;   // LOS vector from receiver antenna to satellite in
                       // local coordinates [km] (n_meas x 3)
  MatrixXd e_rx;       // LOS vector from receiver antenna to satellite in
                       // ECEF coordinates [km] (n_meas x 3)
  VectorXd E;          // Coordinates transformation matrix from the satellite
                       // body‐fixed coordinates to ECEF coordinates
  Vector3d d_rx_disp;  // Displacement by Earth tides at the receiver
                       // position in local coordinates [km] (3)
  VectorXd eps_Phi;    // Phase-range measurement noise [km] (n_bands)

  // Doppler shift measurement
  VectorXd D_rx;  // Doppler shift measurement [Hz] (n_meas * n_bands)
  VectorX r_rx;   // Position of the receiver at time t_rx [km] (3)
  MatrixXd r_tx;  // Position of the transmitter at time t_rx [km] (n_meas x 3)
  MatrixXd v_rx;  // Velocity of the receiver at time t_rx [m/s] (3)
  VectorXd v_tx;  // Velocity of the transmitter at time t_tx [m/s] (n_meas x 3)
  VectorXd eps_D;  // Doppler shift measurement noise [Hz] (n_meas * n_bands)

  // Pseudorange rate measurement
  VectorXd PR_rx;      // Pseudorange rate measurement [m/s] (n_meas * n_bands)
  double dt_rx_dot;    // Receiver clock drift [s/s]
  VectorXd dt_tx_dot;  // Transmitter clock drift [s/s] (n_meas)
  VectorXd eps_PR;     // Pseudorange rate measurement noise [m/s]
                       // (n_meas * n_bands)

  // Link budget
  VectorXd CN0;  // Carrier‐to‐noise density [dB‐Hz] (n_meas * n_bands)

 public:
  GnssMeasurement(const std::vector<Transmission> transmissions);

  GnssMeasurement ExtractSignal(std::string freq_label);

  // Transmission data
  int GetNumMeasurements() const { return n_meas; }
  std::vector<int> GetTxIds() const { return ID_tx; }
  VectorX GetCN0() const { return CN0; }
  VectorXd GetEarthOccultation() const { return vis_earth; }
  VectorXd GetMoonOccultation() const { return vis_moon; }
  VectorXd GetAntennaOccultation() const { return vis_antenna; }
  VectorXd GetIonosOccultation() const { return vis_ionos; }
  VectorXd GetAtmosOccultation() const { return vis_atmos; }

  // Measurement
  VectorX ComputePseudorange(VectorX r_rx, real dt_rx, bool with_noise = false);
  VectorX GetPseudorange(bool with_noise = false);
  VectorX GetPseudorange(double epoch, Vector6 rv_pred, Vector2 clk_pred,
                         MatrixXd &H_pr);
  VectorX GetPseudorange2(double epoch, Vector6 rv_pred, Vector2 clk_pred,
                          MatrixXd &H_pr);
  VectorX GetPseudorangeRate(const VectorX &r_rx_, const VectorX &v_rx_);
  VectorX GetCarrierPhase();
  VectorX GetPhaseRange();
  VectorX GetDopplerShift();

  // Measurement Noise
  double ComputePseudorangeNoise(double CN0);
  double ComputePseudorangeRateNoise(double CN0);
  double ComputeCarrierPhaseNoise(double CN0);
};
}  // namespace lupnt
