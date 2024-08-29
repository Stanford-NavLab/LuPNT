#include <lupnt/physics/spice_interface.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// Demonstrate some basic assertions.
TEST_CASE("SpiceInterface.String2TDB") {
  Real et = spice::String2TDB("2023-04-15 00:00:00 TDB");
  REQUIRE(et.val() == 734788800);
}

TEST_CASE("SpiceInterface.String2TAI") {
  Real t_tai = spice::String2TAI("2023-04-15 00:00:00 TDB");
  REQUIRE(t_tai.val() == 7.347887678143708e+08);
}

TEST_CASE("SpiceInterface.TDBtoStringUTC") {
  Real et = spice::String2TDB("2023-04-15 00:00:00 UTC");
  int prec = 3;
  std::string str = spice::TDBtoStringUTC(et, prec);
  REQUIRE(str == "2023 APR 15 00:00:00.000");
}

TEST_CASE("SpiceInterface.GetBodyPosSpice") {
  Real t_tai = spice::String2TAI("2023-04-15 00:00:00 TDB");
  // 2. GetBodyPosSpice: Get Body Position via SPICE
  auto target = NaifId::MOON;
  auto observer = NaifId::EARTH;
  auto refframe = Frame::GCRF;
  std::string abcorr = "NONE";

  Vec3d pos = spice::GetBodyPosSpice(t_tai, observer, target, refframe, abcorr);
  // std::cout <<  std::setprecision (15) << "Moon Position from Earth at
  // 2023-04-15 00:00:00 TDB (J2000): " << pos[0] << " " << pos[1] << " " <<
  // pos[2] << std::endl;

  double abs_error = 1e-3;
  REQUIRE_THAT(pos[0], WithinAbs(2.636382899441744e+05, abs_error));
  REQUIRE_THAT(pos[1], WithinAbs(-2.210284221463223e+05, abs_error));
  REQUIRE_THAT(pos[2], WithinAbs(-1.318831107682142e+05, abs_error));
}

TEST_CASE("SpiceInterface.GetFrameConversionMat") {
  Real t_tai = spice::String2TAI("2023-04-15 00:00:00 TDB");

  // 3: GetFrameConversionMat
  Mat6 xform(6, 6);
  Mat6 xform_expected{{-0.925133587090201, -0.379635996326334, 0.00208718253372858, 0, 0, 0},
                      {0.379635104286875, -0.925135941018640, -0.000823546386013230, 0, 0, 0},
                      {0.00224357543019373, 3.04773366301969e-05, 0.999997482717042, 0, 0, 0},
                      {2.76834291801792e-05, -6.74619777910846e-05, -5.98953108769266e-08,
                       -0.925133587090201, -0.379635996326334, 0.00208718253372858},
                      {6.74618061429021e-05, 2.76834938727775e-05, -1.52200575820685e-07,
                       0.379635104286875, -0.925135941018640, -0.000823546386013230},
                      {1.47075571708085e-10, 5.94642717972108e-11, -3.31788286899988e-13,
                       0.00224357543019373, 3.04773366301969e-05, 0.999997482717042}};

  xform = spice::GetFrameConversionMat(t_tai, Frame::GCRF, Frame::ITRF);
  RequireNear(xform, xform_expected, 1e-6);
}

TEST_CASE("SpiceInterface.GetBodyPosVel") {
  Real t_tai = spice::String2TAI("2023-04-15 00:00:00 TDB");

  // 4. GetBodyPosVel
  std::string target = "MOON";
  std::string center = "EARTH";

  // NaifId center_id = 399;
  // NaifId target_id = 301;

  Vec6 posvel(6);
  // posvel = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON, Frame::GCRF);
  posvel = spice::GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON);

  Vec6 posvel_expected{263638.289944174,  -221028.422146322, -131883.110768214,
                       0.734154922271287, 0.697461344892098, 0.325673181901724};

  double abs_error = 1e-6;
  RequireNear(posvel, posvel_expected, abs_error);
}
