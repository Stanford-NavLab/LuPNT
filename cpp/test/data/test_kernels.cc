#include <lupnt/conversions/time_conversions.h>
#include <lupnt/core/constants.h>
#include <lupnt/data/kernels.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <iostream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

// Test cases for time conversion
TEST_CASE("data.Kernels") {
  std::filesystem::path test_data_path = GetDataPath().parent_path() / "test" / "data";
  if (!std::filesystem::exists(test_data_path / "body_pos_vel_cspice.txt")) {
    SUCCEED("body_pos_vel_cspice.txt fixture is not available");
    return;
  }

  // Spice
  auto file = OpenTestDataFile("body_pos_vel_cspice.txt");
  std::string line;

  // Time
  std::getline(file, line);
  std::istringstream iss(line);
  std::string tmp;
  double t_tai;
  iss >> tmp >> t_tai;
  Real t_tdb = ConvertTime(Real(t_tai), Time::TAI, Time::TDB);

  // Header
  std::getline(file, line);

  // Values
  while (std::getline(file, line)) {
    std::string center_str, target_str;
    double lt, x, y, z, vx, vy, vz;
    iss.clear();
    iss.str(line);
    iss >> center_str >> target_str >> lt >> x >> y >> z >> vx >> vy >> vz;
    Vec6 rv_in{x, y, z, vx, vy, vz};
    auto center_tmp = enum_cast<BodyId>(center_str);
    auto target_tmp = enum_cast<BodyId>(target_str);
    auto center = enum_cast<BodyId>(center_str).value();
    auto target = enum_cast<BodyId>(target_str).value();
    Vec6 rv = GetBodyPosVel(t_tdb, center, target, Frame::GCRF);
    RequireNear(rv, rv_in, epsilon);

    Vec6 rv_km = GetBodyPosVel(t_tdb, center, target, Frame::GCRF, KM_S_KG_UNITS);
    Vec6 expected_km;
    expected_km.head(3) = rv.head(3) / 1000.0;
    expected_km.tail(3) = rv.tail(3) / 1000.0;
    RequireNear(rv_km, expected_km, epsilon);

    Vec6 rv_tcb_km
        = GetBodyPosVel(t_tdb, center, target, Frame::GCRF, KM_S_KG_UNITS, CoordinateScale::TCB);
    Vec6 expected_tcb_km = expected_km;
    expected_tcb_km.head(3) /= (1.0 - L_B);
    RequireNear(rv_tcb_km, expected_tcb_km, epsilon);
  }
}
