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
#include <lupnt/core/Constants.h>
#include <lupnt/physics/SpiceInterface.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace LPT;
namespace ad = autodiff;
namespace sp = SpiceInterface;

int main() {
  sp::LoadSpiceKernel();

  ad::real et = sp::StringToTDB("2023-04-15 00:00:00 TDB");

  int prec = 3;
  std::string str = sp::TDBtoStringUTC(et, prec);
  std::cout << "TDB: " << et << std::endl;

  ad::real tai = sp::ConvertTime(et, "TDB", "TAI");

  std::cout << "TAI: " << tai << std::endl;

  ad::MatrixXreal xform(6, 6);

  xform = sp::GetFrameConversionMatrix(et, "J2000", "ITRF93");
  std::cout << "XFORM_ITRF: " << std::endl << xform << std::endl;

  xform = sp::GetFrameConversionMatrix(et, "J2000", "IAU_EARTH");
  std::cout << "XFORM_IAUEARTH: " << std::endl << xform << std::endl;

  xform = sp::GetFrameConversionMatrix(et, "J2000", "MOON_PA");
  std::cout << "XFORM_MOONPA: " << std::endl << xform << std::endl;

  xform = sp::GetFrameConversionMatrix(et, "J2000", "IAU_MOON");
  std::cout << "XFORM_IAU: " << std::endl << xform << std::endl;

  return 0;
}