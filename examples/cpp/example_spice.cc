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
namespace ad = autodiff;
namespace sp = SpiceInterface;

int main() {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]");

  sp::LoadSpiceKernel();
  // Print directory from which the executable is being run
  std::cout << "Current working directory: " << std::filesystem::current_path()
            << std::endl;

  real t_tai = sp::StringToTAI("2030-01-01 00:00:00 UTC");
  real t_tdb = sp::ConvertTime(t_tai, TimeSystems::TAI, TimeSystems::TDB);

  int prec = 3;
  std::string str = sp::TAItoStringUTC(t_tai, prec);
  std::cout << "TAI: " << t_tai << std::endl;
  std::cout << "TDB: " << t_tdb << std::endl;

  MatrixX xform(6, 6);
  xform = sp::GetFrameConversionMatrix(t_tai, Frame::GCRF, Frame::ITRF);
  std::cout << "XFORM_ITRF:\n" << xform << std::endl;

  xform = sp::GetFrameConversionMatrix(t_tai, Frame::GCRF, Frame::MOON_PA);
  std::cout << "XFORM_MOONPA:\n" << xform << std::endl;
  xform = sp::GetFrameConversionMatrix(t_tai, Frame::MOON_CI, Frame::MOON_PA);

  Vector3 x =
      sp::GetBodyPos(NaifId::EARTH, t_tai, Frame::GCRF, NaifId::MOON, "NONE");
  std::cout << "r_moon2earth:\n" << x.transpose() << std::endl;

  x = sp::GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH, Frame::GCRF)
          .head(3);
  std::cout << "r_moon2earth:\n" << x.transpose() << std::endl;

  x = sp::GetBodyPos(NaifId::SUN, t_tai, Frame::GCRF, NaifId::MOON, "NONE");
  std::cout << "r_moon2sun:\n" << x.transpose() << std::endl;

  x = sp::GetBodyPosVel(t_tai, NaifId::MOON, NaifId::SUN, Frame::GCRF).head(3);
  std::cout << "r_moon2sun:\n" << x.transpose() << std::endl;

  // Get derivative of r_moon2earth using autodiff
  Vector6 rv_moon2earth;
  auto func = [](const auto &t_tai) {
    return sp::GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH);
  };
  Vector6 rv_moon2earth_dot =
      ad::jacobian(func, wrt(t_tai), at(t_tai), rv_moon2earth);

  std::cout << "v_moon2earth:\n"
            << rv_moon2earth.tail(3).transpose() << std::endl;
  std::cout << "r_moon2earth_dot:\n"
            << rv_moon2earth_dot.head(3).transpose() << std::endl;

  return 0;
}