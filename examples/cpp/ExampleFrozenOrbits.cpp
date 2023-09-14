
#include <lupnt/core/Constants.h>
#include <lupnt/core/UserFilePath.h>
#include <lupnt/dynamics/Dynamics.h>

#include <autodiff/forward/real.hpp>
#include <fstream>
#include <iostream>
#include <lupnt/physics/LunarMeanOsc.cpp>
#include <string>

using namespace LPT;
using namespace std;

typedef ad::VectorXreal Vec;

int main() {
  MoonMeanDynamics moonMeanDynamics;
  ClassicalOE coeMean_1(3738, 0.1, 48.0689 * RAD_PER_DEG, 0, 90 * RAD_PER_DEG,
                        0);
  ClassicalOE coeMean_2(13738, 0.1, 39.6024 * RAD_PER_DEG, 0, 90 * RAD_PER_DEG,
                        0);

  Vec coeMeanVec_1 = coeMean_1.GetVector();
  Vec coeMeanVec_2 = coeMean_2.GetVector();
  Vec delMeanVec_1 = CoeToDelaunay(coeMeanVec_1, MU_MOON, 0, 0);
  Vec delMeanVec_2 = CoeToDelaunay(coeMeanVec_2, MU_MOON, 0, 0);
  Vec coeMeanVec_1_ = DelaunayToCoe(delMeanVec_1, MU_MOON, 0, 0);
  Vec coeMeanVec_2_ = DelaunayToCoe(delMeanVec_2, MU_MOON, 0, 0);

  std::cout << "coeMeanVec_1 : " << coeMeanVec_1.transpose() << std::endl;
  std::cout << "coeMeanVec_1_: " << coeMeanVec_1_.transpose() << std::endl;

  std::cout << "coeMeanVec_2 : " << coeMeanVec_2.transpose() << std::endl;
  std::cout << "coeMeanVec_2_: " << coeMeanVec_2_.transpose() << std::endl;

  Vec coeOscVec_1 = MeanToOsculating(coeMeanVec_1);
  Vec coeOscVec_2 = MeanToOsculating(coeMeanVec_2);

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