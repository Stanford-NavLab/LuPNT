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
#include <lupnt/measurements/occultation.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coord_converter.h>

#include <memory>

using namespace lupnt;

void printOccultation(ad::VectorXreal state_tx_vec,
                      ad::VectorXreal state_rx_vec, std::string seg_planet) {
  ad::VectorXreal tmp_ad(6);
  Eigen::Vector3d segment_eci, segment_mi, user_eci, user_mi;
  double t = 0.0;

  tmp_ad = CoordConverter::Convert(state_tx_vec, t, "GCRF", "MI");
  segment_mi = toEigen(tmp_ad.segment(0, 3));
  tmp_ad = CoordConverter::Convert(state_tx_vec, t, "GCRF", "GCRF");
  segment_eci = toEigen(tmp_ad.segment(0, 3));
  tmp_ad = CoordConverter::Convert(state_rx_vec, t, "MI", "GCRF");
  user_mi = toEigen(tmp_ad.segment(0, 3));
  tmp_ad = CoordConverter::Convert(state_rx_vec, t, "MI", "GCRF");
  user_eci = toEigen(tmp_ad.segment(0, 3));

  std::map<std::string, bool> occ = Occultation::ComputeOccultation(
      segment_eci, segment_mi, user_eci, user_mi, seg_planet);

  for (auto &o : occ) {
    std::cout << o.first << " " << o.second << " ";
  }
  std::cout << std::endl;
}

int main() {
  Eigen::Vector3d segment_eci, segment_mi, user_eci, user_mi;
  ad::VectorXreal state_tx_vec(6), state_rx_vec(6), tmp_ad(6);
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
