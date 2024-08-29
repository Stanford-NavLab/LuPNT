/**
 * @file ExampleOccultation.cpp
 * @author Stanford NAV LAB
 * @brief Example of occultation computation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/occultation.h>

#include <memory>

using namespace lupnt;

void printOccultation(Vec6 state_tx_vec, Vec6 state_rx_vec, std::string seg_planet) {
  Vec6 tmp_ad;
  Vec3d segment_eci, segment_mi, user_eci, user_mi;
  double t = 0.0;

  tmp_ad = ConvertFrame(t, state_tx_vec, Frame::GCRF, Frame::MOON_CI);
  segment_mi = tmp_ad.segment(0, 3).cast<double>();
  tmp_ad = ConvertFrame(t, state_tx_vec, Frame::GCRF, Frame::GCRF);
  segment_eci = tmp_ad.segment(0, 3).cast<double>();
  tmp_ad = ConvertFrame(t, state_rx_vec, Frame::MOON_CI, Frame::GCRF);
  user_mi = tmp_ad.segment(0, 3).cast<double>();
  tmp_ad = ConvertFrame(t, state_rx_vec, Frame::MOON_CI, Frame::GCRF);
  user_eci = tmp_ad.segment(0, 3).cast<double>();

  std::map<std::string, bool> occ
      = Occultation::ComputeOccultationGnss(segment_eci, segment_mi, user_eci, user_mi, seg_planet);

  for (auto &o : occ) {
    std::cout << o.first << " " << o.second << " ";
  }
  std::cout << std::endl;
}

int main() {
  Vec3d segment_eci, segment_mi, user_eci, user_mi;
  Vec6 state_tx_vec, state_rx_vec, tmp_ad;
  double t = 0.0;

  double RE = 6378.0, RM = 1737.0, hGps = 20200.0, hLnss = 5000.0;

  // Satellites
  std::cout << "Earth occult" << std::endl;
  state_tx_vec << -(RE + hGps), 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << -(RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "");

  std::cout << "Moon occult" << std::endl;
  state_tx_vec << (RE + hGps), 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << (RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "");

  std::cout << "Earth moon occult" << std::endl;
  state_tx_vec << -(RE + hGps), 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << (RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "");

  std::cout << "No occult" << std::endl;
  state_tx_vec << (RE + hGps), 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << -(RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "");

  std::cout << "Ionos occult" << std::endl;
  state_tx_vec << (RE + 100.0), 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << -(RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "");

  std::cout << "Tropos occult" << std::endl;
  state_tx_vec << (RE + 10.0), 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << -(RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "");

  // Ground stations
  std::cout << "Earth occult" << std::endl;
  state_tx_vec << 0.0, 0.0, RE, 0.0, 0.0, 0.0;
  state_rx_vec << -(RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "EARTH");

  std::cout << "No occult" << std::endl;
  state_tx_vec << RE, 0.0, 0.0, 0.0, 0.0, 0.0;
  state_rx_vec << -(RM + hLnss), 0.0, 0.0, 0.0, 0.0, 0.0;
  printOccultation(state_tx_vec, state_rx_vec, "EARTH");
}
