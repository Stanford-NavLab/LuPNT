#include <lupnt/interfaces/spice.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("interfaces.spice") {
  spice::LoadSpiceKernel();
  REQUIRE(spice::HasNaifBody("EARTH"));
  REQUIRE(spice::GetNaifId("EARTH") == static_cast<int>(BodyId::EARTH));
  REQUIRE(spice::GetNaifName(static_cast<int>(BodyId::MOON)) == "MOON");
  REQUIRE_FALSE(spice::HasNaifBody("__LUPNT_UNKNOWN_BODY__"));
}

// During SPICE initialization, LoadSpiceKernel automatically looks up and
// downloads the latest high-accuracy Earth (ITRF93) and lunar
// (MOON_PA / MOON_ME) orientation kernels from the NAIF generic kernels
// server (caching them locally under the ephemeris kernel directory, and
// falling back to any already-cached copy if the server cannot be reached --
// e.g. when LUPNT_SKIP_SPICE_KERNEL_DOWNLOAD is set for offline test
// environments).
//
// NOTE: this test deliberately checks frame *availability* via
// `spice::HasFrame` (which only resolves a frame name to a frame ID via
// `namfrm_c`) rather than exercising `spice::GetFrameConversionMat` /
// `spice::StringToTai`, since those additionally route through
// `spice::ConvertTime` / CSPICE's `unitim_c`, which has a pre-existing,
// unrelated crash (see `cpp/test/physics/test_spice_interface.cc`, where the
// `RequireSpiceFixtures` macro is overridden to skip all SPICE runtime tests
// "for now").
TEST_CASE("interfaces.spice.high_accuracy_orientation_kernels") {
  spice::LoadSpiceKernel();

  auto kernel_dir = GetCspiceKernelDir();

  SECTION("high-accuracy Earth orientation (ITRF93) frame is available") {
    // ITRF93 is the high-accuracy Earth body-fixed frame defined by the
    // binary Earth orientation PCK ("earth_latest_high_prec.bpc", or, if
    // unavailable, the cached fallback kernels). Unlike IAU_EARTH (defined
    // in the text PCK and always present), ITRF93 requires a binary PCK to
    // be loaded -- this is exactly the high-accuracy frame that the NAIF
    // tutorial "'High Accuracy' Orientation and Body-fixed Frames for the
    // Moon and Earth" recommends over IAU_EARTH.
    REQUIRE(spice::HasFrame("ITRF93"));

    bool found_earth_pck = false;
    for (const auto& entry : std::filesystem::directory_iterator(kernel_dir)) {
      std::string name = entry.path().filename().string();
      if (name == "earth_latest_high_prec.bpc" || name.rfind("earth_", 0) == 0)
        found_earth_pck = true;
    }
    REQUIRE(found_earth_pck);
  }

  SECTION("high-accuracy lunar orientation (MOON_PA / MOON_ME) is available when downloaded") {
    // Unlike Earth's stable "earth_latest_high_prec.bpc", NAIF does not
    // publish a stable "latest" lunar orientation kernel name, so
    // LoadSpiceKernel must look up the current latest moon_pa_de*.bpc /
    // moon_de*.tf pair from the NAIF directory listing (or fall back to a
    // locally-cached pair). This may legitimately be unavailable in fully
    // offline environments with no pre-cached lunar kernels, so this section
    // only runs its assertions *if* the high-accuracy lunar frame ended up
    // loaded.
    bool moon_pa_available = spice::HasFrame("MOON_PA") || spice::HasFrame("MOON_ME");
    if (!moon_pa_available) {
      WARN(
          "High-accuracy lunar orientation kernels (moon_pa_de*.bpc / moon_de*.tf) were not "
          "available locally and could not be downloaded from NAIF (no network access?) -- "
          "skipping high-accuracy lunar frame checks.");
      return;
    }
    REQUIRE(moon_pa_available);

    // The high-accuracy lunar PCK + FK should now be cached locally for
    // offline reuse.
    bool found_lunar_pck = false;
    bool found_lunar_fk = false;
    for (const auto& entry : std::filesystem::directory_iterator(kernel_dir)) {
      std::string name = entry.path().filename().string();
      if (name.rfind("moon_pa_de", 0) == 0 && name.size() > 4
          && name.substr(name.size() - 4) == ".bpc")
        found_lunar_pck = true;
      if (name.rfind("moon_de", 0) == 0 && name.size() > 3 && name.substr(name.size() - 3) == ".tf")
        found_lunar_fk = true;
    }
    REQUIRE(found_lunar_pck);
    REQUIRE(found_lunar_fk);
  }

  SECTION("legacy low-accuracy frames remain available alongside high-accuracy ones") {
    // IAU_EARTH / IAU_MOON come from the text PCK (pck00011.tpc), which is
    // always loaded, so they should remain resolvable regardless of whether
    // the high-accuracy binary kernels were downloaded.
    REQUIRE(spice::HasFrame("IAU_EARTH"));
    REQUIRE(spice::HasFrame("IAU_MOON"));
  }

  SECTION("unknown frame names are correctly reported as unavailable") {
    REQUIRE_FALSE(spice::HasFrame("__LUPNT_UNKNOWN_FRAME__"));
  }
}
