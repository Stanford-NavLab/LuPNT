/**
 * @file GnssTransmitter.cpp
 * @author Stanford NAV LAB
 * @brief Handles Gnss transmit
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/gnss_transmitter.h"

#include "lupnt/agents/agent.h"
#include "lupnt/core/user_file_path.h"
#include "lupnt/data/kernels.h"
#include "lupnt/measurements/gnss_channel.h"
#include "lupnt/numerics/string_utils.h"

namespace lupnt {

  void GnssTransmitter::InitializeGnssTransmitter() {
    if (gnss_type_ == "GPS") {
      InitializeGPSTransmitter();
    } else if (gnss_type_ == "GLONASS") {
      InitializeGLONASSTransmitter();
    } else if (gnss_type_ == "GALILEO") {
      InitializeGALILEOTransmitter();
    } else if (gnss_type_ == "BEIDOU") {
      InitializeBEIDOUTransmitter();
    } else {
      std::cout << "Gnss type not recognized" << std::endl;
    }
  }

  /**
   * @brief Initialize the GPS transmitter
   *
   */
  void GnssTransmitter::InitializeGPSTransmitter() {
    std::filesystem::path csvpath(GetDataPath() / "gnss" / "gps_table.csv");
    std::vector<std::vector<std::string>> gps_table = ReadCSV(csvpath.string());

    for (size_t i = 0; i < gps_table.size(); i++) {
      if (std::stoi(gps_table[i][0]) == prn_) {
        std::string gps_type = gps_table[i][2];  // 'IIA', 'IIR', 'IIR-M', 'IIF'

        // set transmittion power and antenna pattern, depending on the gps type
        std::string ant_name;
        if (gps_type == "IIA") {
          P_tx = 14.3;                 // dB-W
          ant_name = gps_table[i][4];  // ACE Pattern
          freq_list = {"L1", "L2"};
        } else if (gps_type == "IIR") {
          P_tx = 15.0;                 // dB_W
          ant_name = gps_table[i][3];  // LM Pattern
          freq_list = {"L1", "L2"};
        } else if (gps_type == "IIR-M") {
          P_tx = 15.0;                 // dB_W
          ant_name = gps_table[i][3];  // LM Pattern
          freq_list = {"L1", "L2"};
        } else if (gps_type == "IIF") {
          P_tx = 14.3;                 // dB_W
          ant_name = gps_table[i][4];  // ACE Pattern
          freq_list = {"L1", "L2", "L5"};
        } else if (gps_type == "III") {
          P_tx = 14.3;                 // dB_W
          ant_name = gps_table[i][4];  // ACE Pattern
          freq_list = {"L1", "L2", "L5"};
        } else {
          std::runtime_error("Invalid GPS type");
        }

        std::cout << "PRN: " << prn_ << " type: " << gps_type << " Antenna: " << ant_name
                  << std::endl;

        antenna_ = Antenna(ant_name);
        break;
      }
    }
  }

  /**
   * @brief   Initialize the GLONASS transmitter
   *
   */
  void GnssTransmitter::InitializeGLONASSTransmitter() {
    std::cout << "Antenna type not implemented yet for " << gnss_type_ << std::endl;
  }

  /**
   * @brief Initialize the GALILEO transmitter
   *
   */
  void GnssTransmitter::InitializeGALILEOTransmitter() {
    std::cout << "Antenna type not implemented yet for " << gnss_type_ << std::endl;
  }

  /**
   * @brief Initialize the BEIDOU transmitter
   */
  void GnssTransmitter::InitializeBEIDOUTransmitter() {
    std::cout << "Antenna type not implemented yet for " << gnss_type_ << std::endl;
  }

  /**
   * @brief Get the transmitter orientation
   *        Todo: Change this depending on the gnss type
   *
   * @param t
   * @param r_tx_gcrf
   * @return Vec3d
   */
  std::vector<Vec3d> GnssTransmitter::GetTransmitterOrientation(double t, Vec3d& r_tx_gcrf) {
    // (Sun-Earth) - (Sat-Earth)
    Vec3d r_sat2sun
        = GetBodyPosVel(t, NaifId::EARTH, NaifId::SUN, Frame::GCRF).cast<double>().head(3)
          - r_tx_gcrf;
    auto e_z_gnss = -r_tx_gcrf.normalized();  // Face towards earth center
    auto e_y_gnss = r_sat2sun.cross(r_tx_gcrf).normalized();
    auto e_x_gnss = e_y_gnss.cross(e_z_gnss).normalized();

    std::vector<Vec3d> e_gnss = {e_x_gnss, e_y_gnss, e_z_gnss};
    return e_gnss;
  }

  double GnssTransmitter::GetTransmitterAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) {
    auto e_gnss = GnssTransmitter::GetTransmitterOrientation(t, r_tx_gcrf);
    Vec3d e_x_gnss = e_gnss[0];
    Vec3d e_y_gnss = e_gnss[1];
    Vec3d e_z_gnss = e_gnss[2];
    Vec3d u_tx_rx = (r_rx_gcrf - r_tx_gcrf).normalized();
    double theta_tx = acos(u_tx_rx.dot(e_z_gnss));
    double phi_tx = atan2(u_tx_rx.dot(e_y_gnss), u_tx_rx.dot(e_x_gnss));
    double At = GnssTransmitter::GetAntennaGain(theta_tx * DEG, phi_tx * DEG);
    return At;
  }

  /**
   * @brief Generate a transmission
   *
   * @param t
   * @return Transmission
   */
  GnssTransmission GnssTransmitter::GenerateTransmission(double t) {
    CartesianOrbitState cart_state = GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    ConvertOrbitStateFrame(cart_state, t, Frame::GCRF);

    GnssTransmission trans;
    trans.dt_tx = 0.0;
    trans.r_tx = cart_state.r().cast<double>();
    trans.v_tx = cart_state.v().cast<double>();
    return trans;
  }

}  // namespace lupnt
