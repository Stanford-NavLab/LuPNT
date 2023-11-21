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

#include <autodiff/forward/real.hpp>
#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

typedef Vector6real Vec6;

int main() {
  MoonMeanDynamics moonMeanDynamics;
  ClassicalOE coeMean_1(3738, 0.1, 48.0689 * RAD_PER_DEG, 0, 90 * RAD_PER_DEG,
                        0);
  ClassicalOE coeMean_2(13738, 0.1, 39.6024 * RAD_PER_DEG, 0, 90 * RAD_PER_DEG,
                        0);

  Vec6 coeMeanVec_1 = coeMean_1.GetVector();
  Vec6 coeMeanVec_2 = coeMean_2.GetVector();
  Vec6 delMeanVec_1 = CoeToDelaunay(coeMeanVec_1, MU_MOON, 0, 0);
  Vec6 delMeanVec_2 = CoeToDelaunay(coeMeanVec_2, MU_MOON, 0, 0);
  Vec6 coeMeanVec_1_ = DelaunayToCoe(delMeanVec_1, MU_MOON, 0, 0);
  Vec6 coeMeanVec_2_ = DelaunayToCoe(delMeanVec_2, MU_MOON, 0, 0);

  std::cout << "coeMeanVec_1 : " << coeMeanVec_1.transpose() << std::endl;
  std::cout << "coeMeanVec_1_: " << coeMeanVec_1_.transpose() << std::endl;

  std::cout << "coeMeanVec_2 : " << coeMeanVec_2.transpose() << std::endl;
  std::cout << "coeMeanVec_2_: " << coeMeanVec_2_.transpose() << std::endl;

  Vec6 coeOscVec_1 = MeanToOsculating(coeMeanVec_1, J2_MOON);
  Vec6 coeOscVec_2 = MeanToOsculating(coeMeanVec_2, J2_MOON);

  ClassicalOE coeOsc_1(coeOscVec_1);
  ClassicalOE coeOsc_2(coeOscVec_2);

  coeMean_1.Print();
  cout << endl;
  coeOsc_1.Print();
  cout << endl;
  coeMean_2.Print();
  cout << endl;
  coeOsc_2.Print();
  return 0;
}