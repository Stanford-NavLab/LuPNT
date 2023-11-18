#include <gtest/gtest.h>
#include <lupnt/physics/spice_interface.h>

#include <iostream>

#include "TestUtils.h"

using namespace lupnt;
namespace ad = autodiff;
namespace sp = SpiceInterface;

// Demonstrate some basic assertions.
TEST(SpiceTest, StringToTDB) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");

  // 1.StringToTDB
  ad::real et = sp::StringToTDB("2023-04-15 00:00:00 TDB");
  EXPECT_EQ(et.val(), 734788800);
}

TEST(SpiceTest, StringToTAI) {
  ad::real tai = sp::StringToTAI("2023-04-15 00:00:00 TDB");
  EXPECT_EQ(tai.val(), 7.347887678143708e+08);
}

TEST(SpiceTest, TDBtoStringUTC) {
  ad::real et = sp::StringToTDB("2023-04-15 00:00:00 UTC");
  int prec = 3;
  std::string str = sp::TDBtoStringUTC(et, prec);
  //   std::cout << str << std::endl;
  EXPECT_EQ(str, "2023 APR 15 00:00:00.000");
}

TEST(SpiceTest, GetBodyPos) {
  ad::real et = sp::StringToTDB("2023-04-15 00:00:00 TDB");
  // 2. GetBodyPos: Get Body Position via SPICE
  std::string target = "MOON";
  std::string observer = "EARTH";
  std::string refframe = "J2000";
  std::string abcorr = "NONE";

  Eigen::Vector3d pos = sp::GetBodyPos(target, et, refframe, observer, abcorr);
  // std::cout <<  std::setprecision (15) << "Moon Position from Earth at
  // 2023-04-15 00:00:00 TDB (J2000): " << pos[0] << " " << pos[1] << " " <<
  // pos[2] << std::endl;

  double abs_error = 1e-3;
  EXPECT_NEAR(pos[0], 2.636382899441744e+05, abs_error);
  EXPECT_NEAR(pos[1], -2.210284221463223e+05, abs_error);
  EXPECT_NEAR(pos[2], -1.318831107682142e+05, abs_error);
}

TEST(SpiceTest, GetFrameConversionMatrix) {
  ad::real et = sp::StringToTDB("2023-04-15 00:00:00 TDB");

  // 3: GetFrameConversionMatrix
  ad::MatrixXreal xform(6, 6);
  ad::MatrixXreal xform_expected(6, 6);
  xform_expected << -0.925133587090201, -0.379635996326334, 0.00208718253372858,
      0, 0, 0, 0.379635104286875, -0.925135941018640, -0.000823546386013230, 0,
      0, 0, 0.00224357543019373, 3.04773366301969e-05, 0.999997482717042, 0, 0,
      0, 2.76834291801792e-05, -6.74619777910846e-05, -5.98953108769266e-08,
      -0.925133587090201, -0.379635996326334, 0.00208718253372858,
      6.74618061429021e-05, 2.76834938727775e-05, -1.52200575820685e-07,
      0.379635104286875, -0.925135941018640, -0.000823546386013230,
      1.47075571708085e-10, 5.94642717972108e-11, -3.31788286899988e-13,
      0.00224357543019373, 3.04773366301969e-05, 0.999997482717042;

  xform = sp::GetFrameConversionMatrix(et, "J2000", "ITRF93");
  EXPECT_NEAR_ADMAT(xform, xform_expected, 1e-6);
}

TEST(SpiceTest, GetBodyPosVel) {
  ad::real tai = sp::StringToTAI("2023-04-15 00:00:00 TDB");

  // 4. GetBodyPosVel
  std::string target = "MOON";
  std::string center = "EARTH";

  int center_id = 399;
  int target_id = 301;

  ad::VectorXreal posvel(6);
  posvel = sp::GetBodyPosVel(tai, center_id, target_id);

  ad::VectorXreal posvel_expected(6);
  posvel_expected << 263638.289944174, -221028.422146322, -131883.110768214,
      0.734154922271287, 0.697461344892098, 0.325673181901724;

  double abs_error = 1e-6;
  EXPECT_NEAR_ADVEC(posvel, posvel_expected, abs_error);
}
