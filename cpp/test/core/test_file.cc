#include <lupnt/core/file.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.file") {
  SECTION("CountLines counts newline-delimited records") {
    auto path = std::filesystem::temp_directory_path() / "lupnt_test_count_lines.txt";
    {
      std::ofstream out(path);
      out << "one\n";
      out << "two\n";
      out << "three\n";
    }

    REQUIRE(CountLines(path) == 3);

    std::filesystem::remove(path);
  }

  SECTION("FindFileInDir finds files by filename or stem") {
    auto dir = std::filesystem::temp_directory_path() / "lupnt_test_find_file";
    std::filesystem::create_directories(dir);
    auto path = dir / "target_file.txt";
    {
      std::ofstream out(path);
      out << "content";
    }

    REQUIRE(FindFileInDir(dir, "target_file.txt").has_value());
    REQUIRE(FindFileInDir(dir, "target_file").has_value());
    REQUIRE_FALSE(FindFileInDir(dir, "missing").has_value());

    std::filesystem::remove_all(dir);
  }

  SECTION("PrintDuration formats hours, minutes, and fractional seconds") {
    std::chrono::duration<double> duration(3661.25);

    REQUIRE(PrintDuration(duration, 2) == "1h 1m 1.25s");
  }
}
