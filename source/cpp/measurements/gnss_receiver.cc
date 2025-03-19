/**
 * @file gnss_receiver.cpp
 * @author Stanford NAV LAB
 * @brief Gnss Receiver class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/gnss_receiver.h"

#include "lupnt/data/kernels.h"
#include "lupnt/measurements/gnss_measurement.h"

namespace lupnt {

  /**
   * @brief Initialize the receiver parameters
   *
   */
  void GnssReceiver::InitializeReceiverParams() {
    if (receiver_name_ == "moongpsr") {
      // Lunar Gnss Receiver
      rx_param_.Tsys = 190.0;  // System noise temp [K]
      rx_param_.Ae = 0.0;      // Attenuation due to atmosphere (should be negative) [dB]
      // rx_param_.Nf = -2.85;  // Noise figure of receiver/LNA [dB]
      rx_param_.L = -0.16;            // Receiver implementation, A/D conversion losses [dB]
      rx_param_.As = 0.0;             // System losses, in front of LNA [dB]
      rx_param_.CN0threshold = 15.0;  // CN0 threshold [dB-Hz]

      // GNSS Receiver chip parameter
      gnssr_param_.Bp = 5.0;   // Carrier loop noise bandwidth [Hz]
      gnssr_param_.T = 20e-3;  // Tracking loop integration time [s]
      gnssr_param_.b = 2.0;    // normalized bandwidth [Hz]
      gnssr_param_.Bn = 0.2;   // Code loop noise bandwidth [Hz]
      gnssr_param_.D = 0.3;    // Early-to-late correlator spacing (chips)

    } else {
      std::runtime_error("Receiver name not found");
    }
  }

  /**
   * @brief Get the transmitter orientation
   *        Todo: Change this depending on the gnss type
   *
   * @param t  epoch (TAI)
   * @param r_rx_gcrf  position of the receiver in GCRF [km]
   * @return Vec3d   unit vector of the receiver orientation
   */
  std::vector<Vec3d> GnssReceiver::GetReceiverOrientation(double t, Vec3d& r_rx_gcrf,
                                                          std::string mode) {
    Vec3d r_sat2sun
        = GetBodyPosVel(t, NaifId::EARTH, NaifId::SUN, Frame::GCRF).cast<double>().head(3)
          - r_rx_gcrf;  // (SUN-Earth) - (Sat-Earth) = (Sun-Sat)

    Vec3d e_zero = Vec3d::Zero();

    if (mode == "PZ_EarthPoint") {
      auto e_z = -r_rx_gcrf.normalized();  // Face towards Earth center
      auto e_y = r_sat2sun.cross(r_rx_gcrf).normalized();
      auto e_x = e_y.cross(e_z).normalized();

      std::vector<Vec3d> e_sat = {e_x, e_y, e_z};
      return e_sat;
    } else {
      std::runtime_error("Receiver mode not implemented yet");
      std::vector<Vec3d> e_sat = {e_zero, e_zero, e_zero};
      return e_sat;
    }
  }

  /**
   * @brief Get the Antenna Gain of the Receiver
   *
   * @param t   epoch (TAI) [s]
   * @param r_tx_gcrf  position of the transmitter in GCRF [km]
   * @param r_rx_gcrf  position of the receiver in GCRF [km]
   * @param mode  receiver orientation mode  (PZ_EarthPoint, )
   * @return double
   */
  double GnssReceiver::GetReceiverAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) {
    if (attitude_mode_ == "NONE") {
      std::runtime_error("Receiver attitude mode not set");
    }

    auto e_sat = GnssReceiver::GetReceiverOrientation(t, r_rx_gcrf, attitude_mode_);
    auto e_x = e_sat[0];
    auto e_y = e_sat[1];
    auto e_z = e_sat[2];

    auto u_rx_tx = (r_tx_gcrf - r_rx_gcrf).normalized();
    double theta_rx = acos(u_rx_tx.dot(e_z));
    double phi_rx = atan2(u_rx_tx.dot(e_y), u_rx_tx.dot(e_x));
    double Ar = GnssReceiver::GetAntennaGain(theta_rx * DEG, phi_rx * DEG);
    return Ar;
  }

  /**
   * @brief Get GNSS measurement (without noise) from the receiver at epoch t
   *
   * @param t  receiver epoch (TAI) [s]
   * @return GnssMeasurement
   */
  GnssMeasurement GnssReceiver::GetMeasurement(double t) {
    // Revieve Gnss signals
    std::vector<GnssTransmission> transmissions = gnss_channel_->Receive(*this, t);

    // Generate a measurement from the Gnss transmissions
    if (transmissions.size() == 0) {
      return GnssMeasurement(transmissions);
    }
    GnssMeasurement measurement = GnssMeasurement(transmissions);
    return measurement;
  }

}  // namespace lupnt
