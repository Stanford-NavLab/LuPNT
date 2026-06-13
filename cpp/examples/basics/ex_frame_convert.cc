/**
 * @file ExampleCoordConvert.cpp
 * @author Stanford NAV LAB
 * @brief Compare a Vallado GCRF-to-ITRF frame conversion example against LuPNT.
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <lupnt/conversions/frame_converter.h>
#include <lupnt/core/constants.h>
#include <lupnt/interfaces/spice.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;

int main() {
  auto from = Frame::GCRF;  // J2000
  auto to = Frame::ITRF;    // Earth fixed frame

  // Vallado, p87: reference inertial state at the listed UTC epoch.
  Vec6 rv_gcrf;
  Vec3 pos, vel;
  pos << 5102.5096, 6123.01152, 6378.1368;
  vel << -4.7432196, 0.7905366, 5.553375619;
  rv_gcrf << pos, vel;
  Real t_tdb = spice::StringToTdb("2001/04/06 07:51:28.788 UTC");

  VecX rv_itrf = ConvertFrame(t_tdb, rv_gcrf, from, to);

  std::cout << "rv at J2000 = " << rv_gcrf << std::endl;
  std::cout << "rv at ITRF = " << rv_itrf << std::endl;

  return 0;
}
