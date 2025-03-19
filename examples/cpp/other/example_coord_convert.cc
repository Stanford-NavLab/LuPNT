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
  Vec6 rv_gcrf;
  Vec3 pos, vel;
  pos << 5102.5096, 6123.01152, 6378.1368;
  vel << -4.7432196, 0.7905366, 5.553375619;
  rv_gcrf << pos, vel;
  Real t_tai = spice::String2TAI("2001/04/06 07:51:28.788 UTC");

  VecX rv_itrf = ConvertFrame(t_tai, rv_gcrf, from, to);

  std::cout << "rv at J2000 = " << rv_gcrf << std::endl;
  std::cout << "rv at ITRF = " << rv_itrf << std::endl;

  return 0;
}
