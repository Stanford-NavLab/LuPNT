/**
 * @file ExampleSpice.cpp
 * @author Stanford NAV LAB
 * @brief Example of using the SpiceInterface class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <lupnt/core/constants.h>
#include <lupnt/physics/spice_interface.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace lupnt;

namespace sp = SpiceInterface;

int main() {
  sp::LoadSpiceKernel();
  // Print directory from which the executable is being run
  std::cout << "Current working directory: " << std::filesystem::current_path()
            << std::endl;

  real t_tai = sp::StringToTAI("2022-04-15 00:00:00 TDB");

  int prec = 3;
  std::string str = sp::TAItoStringUTC(t_tai, prec);
  std::cout << "TAI: " << t_tai << std::endl;

  real t_tdb = sp::ConvertTime(t_tai, TimeSystems::TAI, TimeSystems::TDB);

  std::cout << "TDB: " << t_tdb << std::endl;

  MatrixX xform(6, 6);

  xform =
      sp::GetFrameConversionMatrix(t_tai, CoordSystem::GCRF, CoordSystem::ITRF);
  std::cout << "XFORM_ITRF:\n" << xform << std::endl;

  xform =
      sp::GetFrameConversionMatrix(t_tai, CoordSystem::GCRF, CoordSystem::PA);
  std::cout << "XFORM_MOONPA:\n" << xform << std::endl;

  Vector3d x = sp::GetBodyPos(NaifId::EARTH, t_tai, CoordSystem::GCRF,
                              NaifId::MOON, "NONE");
  std::cout << "EARTH2MOON:\n"
            << x.format(Eigen::IOFormat(3, 0, ", ", "\n", "[", "]"))
            << std::endl;
  return 0;
}