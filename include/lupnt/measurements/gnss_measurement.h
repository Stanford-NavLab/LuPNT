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

#include <memory>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/measurements/transmission.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  enum GnssMeasurementType {
    PR,   // Pseudorange
    PRR,  // Pseudorange rate
    CP,   // Carrier phase
  };

  class GnssMeasurement {
    // Implemenation based on Gnss SDR Observables block:
    // https://gnss-sdr.org/docs/sp-blocks/observables/
  private:
    int n_meas = 0;                             // Number of measurements
    std::vector<GnssTransmission> trans_store;  // list of transmittion data
    std::vector<int> ID_tx;                     // ID of the transmitter (n_meas)

    // Visbility
    VecXd vis_earth;
    VecXd vis_moon;
    VecXd vis_antenna;
    VecXd vis_atmos;
    VecXd vis_ionos;

    double c = 299792.458;  // Speed of light [km/s]

    // Pseudorange measurement
    double t_rx;  // Signal reception time measured by the receiver clock [s]
    double t_tx;  // Signal transmission time measured by the transmitter clock [s]

    VecX P_rx;     // Pseudorange measurement [km] (n_meas * n_bands)
    VecXd rho_rx;  // True range from the transmitter to the receiver’s
                   // antenna [km] (n_meas)

    double dt_rx;  // Receiver clock offset from Gnss time [s]
    VecXd dt_tx;   // Transmitter clock offset from Gnss time [s] (n_meas)

    VecXd I_rx;   // Ionospheric delay [km] (n_meas * n_bands)
    VecXd T_rx;   // Tropospheric delay [km] (n_meas)
    VecXd eps_P;  // Pseudorange measurement noise [km] (n_meas)

    // Carrier phase measurement
    VecX phi_rx;      // Carrier phase measurement [cycles] (n_meas * n_bands)
    VecXd phi_rx_tx;  // Phase of the receiver's local oscillator at
                      // time t [cycles] (n_meas * n_bands)
    VecXd phi_tx;     // Phase of the transmitted signal at time t [cycles]
                      // (n_meas * n_bands)
    double t;
    VecXd N_rx;      // Carrier phase int ambiguity [cycles] (n_meas * n_bands)
    VecXd eps_phi;   // Carrier phase measurement noise [cycles]
                     // (n_meas * n_bands)
    VecXd f;         // Carrier frequency [Hz] (n_bands) // TODO: Convert to vector
    double t_0;      // Initial time [s]
    VecXd phi_0;     // Initial phase of the transmitted signal at time t_0
                     // [cycles] (n_meas * n_bands)
    VecXd lambda_;   // Carrier wavelength [km] (n_bands) // TODO: Convert to vector
    VecXd L_rx_txr;  // Integer component of the receiver’s numerically
                     // controlled oscillator (NCO) phase
    VecXd K_rx_txr;  // Integer component of the propagation term

    // Phase-range masurement
    VecXd Phi_rx;     // Phase-range measurement [km] (n_meas * n_bands)
                      // (n_meas * n_bands)
    VecXd B_rx;       // Carrier phase bias [cycles] (n_meas * n_bands)
    VecXd dPhi_rx;    // Receiver’s antenna phase center variation [km]
                      // (n_meas * n_bands)
    VecXd d_rx_pco;   // Receiver’s antenna phase center offset in local
                      // coordinates [km] (n_meas) (n_bands x 3)
    VecXd d_rx_pcv;   // Receiver’s antenna phase center variation [km]
                      // (n_meas * n_bands)
    VecXd d_tx_pco;   // Transmitter’s antenna phase center offset in
                      // local coordinates [km] (n_meas * n_bands x 3)
    VecXd d_tx_pcv;   // Transmitter’s antenna phase center variation
                      // [km] (n_meas * n_bands)
    MatXd e_rx_enu;   // LOS vector from receiver antenna to satellite in
                      // local coordinates [km] (n_meas x 3)
    MatXd e_rx;       // LOS vector from receiver antenna to satellite in
                      // ECEF coordinates [km] (n_meas x 3)
    VecXd E;          // Coordinates transformation matrix from the satellite
                      // body‐fixed coordinates to ECEF coordinates
    Vec3d d_rx_disp;  // Displacement by Earth tides at the receiver
                      // position in local coordinates [km] (3)
    VecXd eps_Phi;    // Phase-range measurement noise [km] (n_bands)

    // Doppler shift measurement
    VecXd D_rx;   // Doppler shift measurement [Hz] (n_meas * n_bands)
    VecX r_rx;    // Position of the receiver at time t_rx [km] (3)
    MatXd r_tx;   // Position of the transmitter at time t_rx [km] (n_meas x 3)
    VecXd v_rx;   // Velocity of the receiver at time t_rx [m/s] (3)
    MatXd v_tx;   // Velocity of the transmitter at time t_tx [m/s] (n_meas x 3)
    VecXd eps_D;  // Doppler shift measurement noise [Hz] (n_meas * n_bands)

    // Pseudorange rate measurement
    VecXd PR_rx;       // Pseudorange rate measurement [m/s] (n_meas * n_bands)
    double dt_rx_dot;  // Receiver clock drift [s/s]
    VecXd dt_tx_dot;   // Transmitter clock drift [s/s] (n_meas)
    VecXd eps_PR;      // Pseudorange rate measurement noise [m/s]
                       // (n_meas * n_bands)

    // Receiver param
    GnssReceiverParam gnssr_param;
    double chip_rate;  // Chip rate [Hz]

    // Link budget
    VecXd CN0;  // Carrier‐to‐noise density [dB‐Hz] (n_meas * n_bands)

  public:
    GnssMeasurement(const std::vector<GnssTransmission> transmissions);

    GnssMeasurement ExtractSignal(std::string freq_label);

    // Transmission data
    int GetTrackedSatelliteNum() const { return n_meas; }

    std::vector<int> GetTxIds() const { return ID_tx; }
    VecX GetCN0() const { return CN0; }
    VecXd GetEarthOccultation() const { return vis_earth; }
    VecXd GetMoonOccultation() const { return vis_moon; }
    VecXd GetAntennaOccultation() const { return vis_antenna; }
    VecXd GetIonosOccultation() const { return vis_ionos; }
    VecXd GetAtmosOccultation() const { return vis_atmos; }

    /***********************************************************
     * General Methods for computing Measurements
     ***********************************************************/

    /**
     * @brief Compute the pseudorange measurement
     *
     * @param r_rx   Receiver position [km]
     * @param dt_rx  Receiver clock offset from Gnss time [s]
     * @param with_noise   use noise
     * @param seed   random seed
     * @return VecX
     */
    VecX ComputeGnssPseudorange(VecX r_rx, Real dt_rx, bool with_noise = false, int seed = 0);
    VecX ComputeGnssPseudorangerate(VecX r_rx, VecX v_rx, Real dt_rx_dot, bool with_noise = false,
                                    int seed = 0);
    VecX ComputeGnssCarrierPhase(VecX r_rx, Real dt_rx, VecX N_rx, bool with_noise = false,
                                 int seed = 0);

    /***********************************************************
     *  Methods for true measurement generation
     ***********************************************************/

    /**
     * @brief Get the Gnss Measurement for the observed signal
     *
     * @param meas_type  vector of measurement types
     * @param with_noise  use noise
     * @param seed   random seed
     * @return VecX
     */
    VecX GetGnssMeasurement(std::vector<GnssMeasurementType> meas_type, bool with_noise = false,
                            int seed = 0);
    /**
     * @brief Get the Pseudorange for the observed signal
     *
     * @param with_noise  use noise
     * @param seed  random seed
     * @return VecX
     */
    VecX GetPseudorange(bool with_noise = false, int seed = 0);

    /**
     * @brief Get the Pseudorange rate for the observed signal
     *
     * @param with_noise  use noise
     * @param seed  random seed
     * @return VecX
     */
    VecX GetPseudorangerate(bool with_noise = false, int seed = 0);

    /**
     * @brief Get the Carrier Phase for the observed signal
     *
     * @param with_noise  use noise
     * @param seed  random seed
     * @return VecX
     */
    VecX GetCarrierPhase(bool with_noise = false, int seed = 0);

    /***********************************************************
     *  Methods for predicted measurement generation
     ***********************************************************/
    /**
     * @brief Get the Gnss Measurement for the predicted state
     *       This method is used for the measurement update step
     *       in the filter
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param frame_in  coordinate system of the input state
     */
    VecX GetPredictedGnssMeasurement(double epoch, Vec6 rv_pred, Vec2 clk_pred, VecX N_pred,
                                     MatXd &H_gnss, std::vector<GnssMeasurementType> meas_type,
                                     Frame frame_in = Frame::MOON_CI);

    /**
     * @brief Get the Pseudorange for the predicted state
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_pr       Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return VecX
     */
    VecX GetPredictedPseudorange(double epoch, Vec6 rv_pred, Vec2 clk_pred, MatXd &H_pr,
                                 Frame frame_in = Frame::MOON_CI);

    /**
     * @brief Get the Pseudorange Analytical Jacobian object
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_pr      Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return * VecX
     */
    VecX GetPredictedPseudorangeAnalyticalJacobian(double epoch, Vec6 rv_pred, Vec2 clk_pred,
                                                   MatXd &H_pr, Frame frame_in = Frame::MOON_CI);

    /**
     * @brief Get the Pseudorange Rate object
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_prr     Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return VecX
     */
    VecX GetPredictedPseudorangerate(double epoch, Vec6 rv_pred, Vec2 clk_pred, MatXd &H_prr,
                                     Frame frame_in = Frame::MOON_CI);

    /**
     * @brief Get the Carrier Phase object
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_cp      Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return VecX
     */
    VecX GetPredictedCarrierPhase(double epoch, Vec6 rv_pred, Vec2 clk_pred, VecX N_pred,
                                  MatXd &H_cp, Frame frame_in = Frame::MOON_CI);

    /*********************************************************************
     * Noise Models
     ********************************************************************/
    VecXd GetGnssNoiseStdVec(std::vector<GnssMeasurementType> meas_type);
    VecXd GetPseudorangeNoiseStdVec();
    VecXd GetPseudorangeRateNoiseStdVec();
    VecXd GetCarrierPhaseNoiseStdVec();

    void SetGnssReceiverParam(GnssReceiverParam gnssr_param_input) {
      gnssr_param = gnssr_param_input;
    }

    /**
     * @brief Compute the pseudorange noise using thermal noise in DLL
     * Reference: Reference: "Understanding GPS", p195
     *
     * @param CN0  Carrier‐to‐noise density [dB‐Hz]
     * @return double  Pseudorange noise [km]
     */
    double ComputeGnssPseudorangeNoise(double CN0);

    /**
     * @brief Compute the pseudorange rate noise using thermal noise in FLL
     * Reference: "Understanding GPS", p192
     *
     * @param CN0  Carrier‐to‐noise density [dB‐Hz]
     * @return double  Pseudorange rate noise [km/s]
     */
    double ComputeGnssPseudorangerateNoise(double CN0, double lambda);

    /**
     * @brief Compute the carrier phase noise using thermal noise in PLL
     * Reference: "Understanding GPS", p185
     *
     * @param CN0  Carrier‐to‐noise density [dB‐Hz]
     * @return double  Carrier phase noise [cycles]
     */
    double ComputeGnssCarrierPhaseNoise(double CN0, double lambda);
  };
}  // namespace lupnt
