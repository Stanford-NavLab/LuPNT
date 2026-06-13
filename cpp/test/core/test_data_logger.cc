#include <lupnt/core/data_logger.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <highfive/H5Easy.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.data_logger") {
  std::filesystem::path path = std::filesystem::temp_directory_path() / "lupnt_data_logger_test.h5";

  DataLogger::SetOutputFile(path);
  DataLogger::Log("data_logger_test/scalar", 1.25);
  DataLogger::Log("data_logger_test/real", Real(2.5));
  DataLogger::Log("data_logger_test/string", std::string("ok"));

  Vec3d vec(1.0, 2.0, 3.0);
  DataLogger::Log("data_logger_test/vector", vec);

  ClockState3 state(Vec3(4.0, 5.0, 6.0));
  DataLogger::LogState("data_logger_test/state", state);
  DataLogger::Flush();

  H5Easy::File file(path, H5Easy::File::ReadOnly);
  auto scalar = H5Easy::load<std::vector<double>>(file, "data_logger_test/scalar");
  auto real = H5Easy::load<std::vector<double>>(file, "data_logger_test/real");
  auto text = H5Easy::load<std::vector<std::string>>(file, "data_logger_test/string");
  auto matrix = H5Easy::load<std::vector<MatXd>>(file, "data_logger_test/vector");
  auto bias = H5Easy::load<std::vector<double>>(file, "data_logger_test/state/b_s");

  REQUIRE(scalar.size() == 1);
  REQUIRE(real.size() == 1);
  REQUIRE(text.size() == 1);
  REQUIRE(matrix.size() == 1);
  REQUIRE(bias.size() == 1);
  REQUIRE_THAT(scalar.front(), Catch::Matchers::WithinAbs(1.25, epsilon));
  REQUIRE_THAT(real.front(), Catch::Matchers::WithinAbs(2.5, epsilon));
  REQUIRE(text.front() == "ok");
  REQUIRE(matrix.front().isApprox(vec, epsilon));
  REQUIRE_THAT(bias.front(), Catch::Matchers::WithinAbs(4.0, epsilon));
}
