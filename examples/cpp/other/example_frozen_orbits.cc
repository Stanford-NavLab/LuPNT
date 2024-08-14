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
#include <lupnt/physics/orbit_state.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

int main() {
  MoonMeanDynamics moonMeanDynamics;

  Vec6 coe_m_1{3738, 0.1, 48.0689 * RAD, 0, 90 * RAD, 0};
  Vec6 coe_m_2{13738, 0.1, 39.6024 * RAD, 0, 90 * RAD, 0};

  auto deloe_m_1 = Classical2Delaunay(coe_m_1, GM_MOON);
  auto deloe_m_2 = Classical2Delaunay(coe_m_2, GM_MOON);
  auto coe_m_1_ = Delaunay2Classical(deloe_m_1, GM_MOON);
  auto coe_m_2_ = Delaunay2Classical(deloe_m_2, GM_MOON);

  std::cout << "coe_m_1 : " << coe_m_1.transpose() << std::endl;
  std::cout << "coe_m_1_: " << coe_m_1_.transpose() << std::endl;

  std::cout << "coe_m_2 : " << coe_m_2.transpose() << std::endl;
  std::cout << "coe_m_2_: " << coe_m_2_.transpose() << std::endl;

  auto coe_c_1 = Mean2Osculating(coe_m_1, GM_MOON, J2_MOON);
  auto coe_c_2 = Mean2Osculating(coe_m_2, GM_MOON, J2_MOON);
  return 0;
}
