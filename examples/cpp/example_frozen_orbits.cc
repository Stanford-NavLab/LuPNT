/**
 * @file ExampleFrozenOrbits.cpp
 * @author Stanford NAV LAB
 * @brief Example of frozen orbits
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <lupnt/core/constants.h>
#include <lupnt/core/user_file_path.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/physics/orbit_state_utils.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

int main() {
  MoonMeanDynamics moonMeanDynamics;

  Vector6real coe_m_1{3738, 0.1, 48.0689 * RAD_PER_DEG, 0, 90 * RAD_PER_DEG, 0};
  Vector6real coe_m_2{13738, 0.1, 39.6024 * RAD_PER_DEG, 0, 90 * RAD_PER_DEG,
                      0};

  auto deloe_m_1 = CoeToDelaunay(coe_m_1, MU_MOON, 0, 0);
  auto deloe_m_2 = CoeToDelaunay(coe_m_2, MU_MOON, 0, 0);
  auto coe_m_1_ = DelaunayToCoe(deloe_m_1, MU_MOON, 0, 0);
  auto coe_m_2_ = DelaunayToCoe(deloe_m_2, MU_MOON, 0, 0);

  std::cout << "coe_m_1 : " << coe_m_1.transpose() << std::endl;
  std::cout << "coe_m_1_: " << coe_m_1_.transpose() << std::endl;

  std::cout << "coe_m_2 : " << coe_m_2.transpose() << std::endl;
  std::cout << "coe_m_2_: " << coe_m_2_.transpose() << std::endl;

  auto coe_c_1 = MeanToOsc(coe_m_1, J2_MOON);
  auto coe_c_2 = MeanToOsc(coe_m_2, J2_MOON);
  return 0;
}