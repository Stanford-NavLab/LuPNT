#include <lupnt/conversions/attitude_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.attitude_conversions") {
  SECTION("scalar-first and scalar-last quaternion conventions round trip") {
    Quaternion q_sf(Vec4(0.5, 0.5, 0.5, 0.5));

    State q_sl = ScalarFirstToLast(q_sf);
    State recovered = ScalarLastToFirst(q_sl);

    REQUIRE_THAT(q_sl(0).val(), WithinAbs(0.5, epsilon));
    REQUIRE_THAT(q_sl(3).val(), WithinAbs(0.5, epsilon));
    for (int i = 0; i < 4; ++i) REQUIRE_THAT(recovered(i).val(), WithinAbs(q_sf(i).val(), epsilon));
  }

  SECTION("quaternion and rotation matrix conversions round trip") {
    Mat3 R = RollPitchYawToRot(Vec3(0.1, -0.2, 0.3));
    State q = RotToQuat(R);
    Mat3 recovered = QuatToRot(q);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        REQUIRE_THAT(recovered(i, j).val(), WithinAbs(R(i, j).val(), 1.0e-12));
  }

  SECTION("roll pitch yaw conversions round trip") {
    Vec3 rpy(0.3, -0.2, 0.1);

    Vec3 recovered = RotToRollPitchYaw(RollPitchYawToRot(rpy));

    for (int i = 0; i < 3; ++i) REQUIRE_THAT(recovered(i).val(), WithinAbs(rpy(i).val(), epsilon));
  }
}
