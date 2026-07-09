#include <lupnt/devices/gnss_device.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("devices.space_comms") {
  SECTION("GNSS frequency and chip-rate maps expose common signals") {
    REQUIRE_THAT(GNSS_FREQ_MAP.at(GnssFreq::L1).val(), WithinAbs(1575.42e6, 1.0));
    REQUIRE_THAT(GNSS_FREQ_MAP.at(GnssFreq::L5).val(), WithinAbs(1176.45e6, 1.0));
    REQUIRE_THAT(GNSS_RC_MAP.at(GnssFreq::L1).val(), WithinAbs(1.023e6, 1.0));
    REQUIRE_THAT(GNSS_RC_MAP.at(GnssFreq::E5a).val(), WithinAbs(10.23e6, 1.0));
  }

  SECTION("Galileo transmitter initializes all supported frequencies") {
    GnssTransmitter tx(GnssConst::GALILEO, 1);
    auto freqs = tx.GetFreqList();

    REQUIRE(freqs.size() == 4);
    REQUIRE(freqs[0] == GnssFreq::E1);
    REQUIRE(freqs[1] == GnssFreq::E5a);
    REQUIRE(freqs[2] == GnssFreq::E5b);
    REQUIRE(freqs[3] == GnssFreq::E6);
    for (auto freq : freqs) {
      REQUIRE(tx.GetTransmitPower().count(freq) == 1);
      REQUIRE(tx.GetChipRate().count(freq) == 0);
      REQUIRE(tx.GetAntennas().count(freq) == 1);
    }
  }

  SECTION("GPS transmitter powers follow block_power reference mapping") {
    struct Case {
      int prn;
      double l1_dbw;
      double l5_dbw;
    };

    const std::vector<Case> cases = {
        {13, 17.3, -82.7},  // IIR
        {5, 18.8, -81.2},   // IIR_M
        {1, 16.2, 19.2},    // IIF
        {4, 18.8, 21.8},    // III
    };

    for (const auto& c : cases) {
      GnssTransmitter tx(GnssConst::GPS, c.prn);
      auto powers = tx.GetTransmitPower();
      REQUIRE_THAT(powers.at(GnssFreq::L1).val(), WithinAbs(c.l1_dbw, epsilon));
      REQUIRE_THAT(powers.at(GnssFreq::L5).val(), WithinAbs(c.l5_dbw, epsilon));
    }
  }
}
