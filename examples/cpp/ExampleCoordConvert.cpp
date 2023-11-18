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
#include <lupnt/physics/coord_converter.h>
#include <lupnt/physics/spice_interface.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;
namespace ad = autodiff;
namespace sp = SpiceInterface;

int main() {
  double t_eph = 0.0;  // seconds since J2000

  // Get the position and velocity of the target body relative to the center
  // body
  std::string from = "GCRF";  // J2000
  std::string to = "ITRF";    // Earth fixed frame

  // Vallado, p87
  ad::VectorXreal posvel_GCRF(6), pos(3), vel(3);
  pos << 5102.5096, 6123.01152, 6378.1368;
  vel << -4.7432196, 0.7905366, 5.553375619;
  posvel_GCRF << pos, vel;
  ad::real tai = sp::StringToTAI("2001/04/06 07:51:28.788 UTC");

  ad::VectorXreal posvel_ITRF =
      CoordConverter::Convert(posvel_GCRF, tai, from, to);

  std::cout << "Posvel at J2000 = " << posvel_GCRF << std::endl;
  std::cout << "Posvel at ITRF = " << posvel_ITRF << std::endl;

  return 0;
}