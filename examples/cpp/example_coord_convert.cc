/**
 * @file ExampleCoordConvert.cpp
 * @author Stanford NAV LAB
 * @brief Example for coordinate conversion
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <lupnt/core/constants.h>
#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/spice_interface.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;

int main() {
  double t_eph = 0.0;  // seconds since J2000

  // Get the position and velocity of the target body relative to the center
  // body
  auto from = Frame::GCRF;  // J2000
  auto to = Frame::ITRF;    // Earth fixed frame

  // Vallado, p87
  Vec6 posvel_GCRF;
  Vec3 pos, vel;
  pos << 5102.5096, 6123.01152, 6378.1368;
  vel << -4.7432196, 0.7905366, 5.553375619;
  posvel_GCRF << pos, vel;
  real tai = StringToTAI("2001/04/06 07:51:28.788 UTC");

  VecX posvel_ITRF = FrameConverter::Convert(tai, posvel_GCRF, from, to);

  std::cout << "Posvel at J2000 = " << posvel_GCRF << std::endl;
  std::cout << "Posvel at ITRF = " << posvel_ITRF << std::endl;

  return 0;
}