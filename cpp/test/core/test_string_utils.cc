#include <lupnt/core/string_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.string_utils") {
  SECTION("SplitString preserves empty fields") {
    auto fields = SplitString("alpha,beta,,delta", ',');

    REQUIRE(fields == std::vector<std::string>{"alpha", "beta", "", "delta"});
  }

  SECTION("ReadCsv skips the header and trims whitespace") {
    auto path = std::filesystem::temp_directory_path() / "lupnt_test_string_utils.csv";
    {
      std::ofstream out(path);
      out << "name, value\n";
      out << " first , 1 \n";
      out << "second,2\r\n";
    }

    auto rows = ReadCsv(path);

    REQUIRE(rows.size() == 2);
    REQUIRE(rows[0] == std::vector<std::string>{"first", "1"});
    REQUIRE(rows[1] == std::vector<std::string>{"second", "2"});

    std::filesystem::remove(path);
  }
}
