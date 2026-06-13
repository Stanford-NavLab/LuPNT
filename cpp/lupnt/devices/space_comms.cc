#include "lupnt/devices/space_comms.h"

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/core/string_utils.h"

namespace lupnt {

  bool StringInCandidates(std::vector<std::string>& candidates, std::string str) {
    for (size_t i = 0; i < candidates.size(); i++)
      if (str.find(candidates[i]) != std::string::npos) return true;
    return false;
  }

  const std::map<GnssFreq, Real> GNSS_FREQ_MAP
      = {{GnssFreq::L1, 1575.42e6},  {GnssFreq::L2, 1227.60e6}, {GnssFreq::L5, 1176.45e6},
         {GnssFreq::E1, 1575.42e6},  {GnssFreq::E6, 1278.75e6}, {GnssFreq::E5, 1191.795e6},
         {GnssFreq::E5a, 1176.45e6}, {GnssFreq::E5b, 1207.14e6}};

  const std::map<GnssFreq, Real> GNSS_RC_MAP
      = {{GnssFreq::L1, 1.023e6},  {GnssFreq::L2, 0.5115e6}, {GnssFreq::L5, 10.23e6},
         {GnssFreq::E1, 1.023e6},  {GnssFreq::E6, 0.5115e6}, {GnssFreq::E5, 10.23e6},
         {GnssFreq::E5a, 10.23e6}, {GnssFreq::E5b, 10.23e6}};

  namespace {
    Real GpsBlockTransmitPowerDbw(const std::string& block_name, GnssFreq freq) {
      Real freq_offset = 0.0;
      if (freq == GnssFreq::L5) {
        if (block_name == "IIF" || block_name == "III") {
          freq_offset = 3.0;
        } else if (block_name == "IIR" || block_name == "IIR_M") {
          freq_offset = -100.0;
        }
      }

      if (block_name == "IIR") return 17.3 + freq_offset;
      if (block_name == "IIR_M") return 18.8 + freq_offset;
      if (block_name == "IIF") return 16.2 + freq_offset;
      if (block_name == "III") return 18.8 + freq_offset;
      return 14.0 + freq_offset;
    }
  }  // namespace

  // GnssTransmitter **************************************************************

  // Constructor
  GnssTransmitter::GnssTransmitter(GnssConst gnss_const, int prn)
      : Transmitter(), gnss_const_(gnss_const), prn_(prn) {
    switch (gnss_const_) {
      case GnssConst::GPS: InitGps(); break;
      case GnssConst::GLONASS: InitGlonass(); break;
      case GnssConst::GALILEO: InitGalileo(); break;
      case GnssConst::BEIDOU: InitBeidou(); break;
      case GnssConst::QZSS: InitQzss(); break;
      default: LUPNT_CHECK(false, "Invalid GNSS type", "GnssTransmitter");
    }
  }

  // Constructor from YAML
  GnssTransmitter::GnssTransmitter(Config& config) : Transmitter(config) {
    gnss_const_ = enum_cast<GnssConst>(config["gnss_const"].as<std::string>()).value();
    prn_ = config["prn"].as<int>();

    switch (gnss_const_) {
      case GnssConst::GPS: InitGps(); break;
      case GnssConst::GLONASS: InitGlonass(); break;
      case GnssConst::GALILEO: InitGalileo(); break;
      case GnssConst::BEIDOU: InitBeidou(); break;
      case GnssConst::QZSS: InitQzss(); break;
      default: LUPNT_CHECK(false, "Invalid GNSS type", "GnssTransmitter");
    }
  }

  // GPS
  void GnssTransmitter::InitGps() {
    // Todo: Autonomously generate the prn to svn mapping
    std::filesystem::path csvpath = GetFilePath("gps_table.csv");
    std::vector<std::vector<std::string>> gps_table = ReadCsv(csvpath);

    for (size_t i = 0; i < gps_table.size(); i++) {
      if (std::stoi(gps_table[i][0]) != prn_) continue;

      // set transmittion power and antenna pattern, depending on the gps type
      auto gps_type = gps_table[i][2];

      std::string ant_name;
      if (gps_type == "IIA") {
        ant_name = gps_table[i][4];  // ACE Pattern
        freq_list_ = {GnssFreq::L1, GnssFreq::L2};

      } else if (gps_type == "IIR" || gps_type == "IIR_M") {
        ant_name = gps_table[i][3];  // LM Pattern
        freq_list_ = {GnssFreq::L1, GnssFreq::L5};
      } else if (gps_type == "IIF") {
        ant_name = gps_table[i][4];  // ACE Pattern
        freq_list_ = {GnssFreq::L1, GnssFreq::L5};
      } else if (gps_type == "III") {
        ant_name = gps_table[i][3];  // LM Pattern
        freq_list_ = {GnssFreq::L1, GnssFreq::L5};
      }

      for (auto freq : freq_list_) P_tx_[freq] = GpsBlockTransmitPowerDbw(gps_type, freq);

      if (gps_type == "III") {
        // For Block III, the antenna pattern is different for each frequency
        for (auto freq : freq_list_)
          antennas_[freq] = Antenna(ant_name + "_" + std::string(enum_name(freq)) + ".txt");
      } else {
        // For other GPS types, we assume the antenna pattern is the same for all frequencies
        for (auto freq : freq_list_) antennas_[freq] = Antenna(ant_name);
      }

      return;
    }
    LUPNT_CHECK(false, "PRN not found in GPS table", "GnssTransmitter");
  }

  // GLONASS
  void GnssTransmitter::InitGlonass() {
    LUPNT_CHECK(false, "GLONASS not implemented yet", "GnssTransmitter");
  }

  // Galileo
  void GnssTransmitter::InitGalileo() {
    freq_list_ = {GnssFreq::E1, GnssFreq::E5a, GnssFreq::E5b, GnssFreq::E6};
    for (auto freq : freq_list_) {
      antennas_[freq] = Antenna("Galileo_" + std::string(enum_name(freq)) + ".txt");
    }
    for (auto freq : freq_list_) {
      P_tx_[freq] = 14.0;  // [dB-W]
    }
  }

  // Beidou
  void GnssTransmitter::InitBeidou() {
    LUPNT_CHECK(false, "BEIDOU not implemented yet", "GnssTransmitter");
  }

  // QZSS
  void GnssTransmitter::InitQzss() {
    std::vector<std::string> names = {"1R", "02", "03", "04", "05", "06", "07"};
    if (prn_ <= 4) {
      freq_list_ = {GnssFreq::L1, GnssFreq::L2, GnssFreq::L5};
    } else {
      freq_list_ = {GnssFreq::L1, GnssFreq::L5};
    }
    for (auto freq : freq_list_) {
      antennas_[freq]
          = Antenna("QZSS_" + names[prn_ - 1] + "_" + std::string(enum_name(freq)) + ".txt");
    }
    // Reference: Enhancing Navigation Accuracy in a Geostationary Orbit by Utilizing a Regional
    // Navigation Satellite System
    for (auto freq : freq_list_) {
      P_tx_[freq] = 14.1;  // [dB-W]
    }
  }

  REGISTER_FACTORY_CLASS(Device, GnssTransmitter)

  // GnssReceiver **************************************************************

  // Constructor
  GnssReceiver::GnssReceiver(Config& config) : Receiver(config) {}

  REGISTER_FACTORY_CLASS(Device, GnssReceiver)

}  // namespace lupnt
