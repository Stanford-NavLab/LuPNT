// Regression tests for two "declared but broken" API fixes:
//   * filter_utils GetProcessNoiseFunction was declared (identity pass-through)
//     with no definition -> unlinkable for any caller. Now defined.
//   * string_file_utils get_file_path was declared with `const std::string&` in
//     the header but defined with `std::string_view` in the .cc, so the declared
//     overload had no definition (a latent link error). The header now matches.

#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "../utils.cc"
#include "lupnt/environment/plasma/core/string_file_utils.h"
#include "lupnt/numerics/filters/filter_utils.h"

using namespace lupnt;

TEST_CASE("numerics.filter_utils.process_noise_identity") {
  ProcessNoiseFunction q
      = [](const State&, Real, Real) -> MatXd { return MatXd::Identity(3, 3) * 2.0; };
  ProcessNoiseFunction wrapped = GetProcessNoiseFunction(q);
  MatXd out = wrapped(State(), Real(0.0), Real(1.0));
  REQUIRE(out.rows() == 3);
  REQUIRE(out.cols() == 3);
  REQUIRE(out(1, 1) == 2.0);  // identity pass-through returns the same function
}

TEST_CASE("environment.plasma.get_file_path_resolves_bundled") {
  pecsim::set_base_path((GetDataPath() / "plasma").string());
  // Calling get_file_path via the header must link (the fix) and resolve a
  // bundled solar-index file under <base>/data.
  std::filesystem::path p = pecsim::get_file_path("apf107.dat");
  REQUIRE(p.filename() == "apf107.dat");
  REQUIRE(std::filesystem::exists(p));
}
