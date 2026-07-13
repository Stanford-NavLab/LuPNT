#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <vector>

#include "lupnt/core/file.h"
#include "lupnt/environment/plasma/core/constants.h"
#include "lupnt/environment/plasma/core/user_filepath.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"
#include "lupnt/environment/plasma/gcpm/gcpm_interface.h"
#include "lupnt/environment/plasma/gcpm/iri_interface.h"
#include "lupnt/environment/plasma/igrf/igrf_interface.h"

using namespace pecsim;
using namespace Catch::Matchers;

// The IRI/GCPM code loads coefficient files relative to a base path and calls
// std::filesystem::current_path() internally. We point the base path at the
// bundled plasma data dir and restore the working directory afterwards so no
// other test in the binary is affected.
static void SetupPlasmaBasePath() {
  std::filesystem::path plasma_base = lupnt::GetDataPath() / "plasma";
  set_base_path(plasma_base.string());
}

namespace {
  // Restores the process working directory on scope exit (the GCPM/IRI entry
  // points chdir into the data directory while running).
  struct CwdGuard {
    std::filesystem::path saved;
    CwdGuard() : saved(std::filesystem::current_path()) {}
    ~CwdGuard() { std::filesystem::current_path(saved); }
  };

  // A representative daytime, mid-solar-cycle epoch that is inside the range of
  // the bundled apf107 / ig_rz data files.
  DateTime SampleEpoch() { return DateTime{2015, 172, 12, 0, 0.0}; }

  double Magnitude3(const std::vector<double>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }
}  // namespace

// -------------------------------- constants_gcpm ---------------------------

TEST_CASE("environment.plasma.gcpm.constants_pn_ps") {
  // PN / PS index into the auroral-edge tables; the first entry is a documented
  // constant (origPN[0] == 73.6) and every entry is an invariant latitude in a
  // physically sane band.
  REQUIRE_THAT(PN(0, 0), WithinAbs(73.6, 1e-9));
  REQUIRE(std::isfinite(PS(0, 0)));
  for (int j = 0; j < 10; ++j) {
    for (int i = 0; i < 72; ++i) {
      REQUIRE(std::abs(PN(i, j)) <= 90.0);
      REQUIRE(std::abs(PS(i, j)) <= 90.0);
    }
  }
}

// ----------------------------------- IGRF ----------------------------------

TEST_CASE("environment.plasma.gcpm.igrf14_field") {
  CwdGuard guard;
  SetupPlasmaBasePath();

  const double lat_rad = 0.35;  // ~20 deg N
  const double lon_rad = 0.5;   // ~29 deg E
  const double year = 2015.5;
  const double r_surface_km = 6371.2;

  std::vector<double> B_surface = igrf14(lat_rad, lon_rad, r_surface_km, year);
  REQUIRE(B_surface.size() == 3);
  for (double c : B_surface) REQUIRE(std::isfinite(c));

  // Total field near the surface is within a broad, physically-motivated band
  // (roughly 20,000-70,000 nT anywhere on Earth).
  double f_surface = Magnitude3(B_surface);
  REQUIRE(f_surface > 1.0e4);
  REQUIRE(f_surface < 1.0e5);

  // The field weakens with geocentric distance (dipole ~ 1/r^3), so a point
  // 2000 km higher must have a strictly smaller magnitude.
  std::vector<double> B_high = igrf14(lat_rad, lon_rad, r_surface_km + 2000.0, year);
  double f_high = Magnitude3(B_high);
  REQUIRE(std::isfinite(f_high));
  REQUIRE(f_high < f_surface);
}

// ----------------------------------- IRI -----------------------------------

TEST_CASE("environment.plasma.gcpm.iri_2020_density") {
  CwdGuard guard;
  SetupPlasmaBasePath();

  DateTime dt = SampleEpoch();
  const double amlt = 12.0;  // local noon
  const double alatr = 0.2;  // ~11 deg magnetic latitude
  const double akp = 3.0;

  // ~400 km altitude (r in Earth radii); solidly in the F-region.
  double r_low = 1.0 + 400.0 / RE;
  std::vector<double> n_low = iri_2020(dt, r_low, amlt, alatr, akp);
  REQUIRE(n_low.size() == 4);
  for (double c : n_low) {
    REQUIRE(std::isfinite(c));
    REQUIRE(c >= 0.0);
  }
  // Daytime F-region electron density is positive.
  REQUIRE(n_low[0] > 0.0);

  // ~1500 km is well into the topside, where electron density decreases with
  // altitude, so the higher point has no more electrons than the lower one.
  double r_high = 1.0 + 1500.0 / RE;
  std::vector<double> n_high = iri_2020(dt, r_high, amlt, alatr, akp);
  REQUIRE(n_high[0] >= 0.0);
  REQUIRE(n_low[0] >= n_high[0]);
}

TEST_CASE("environment.plasma.gcpm.iri_2007_density") {
  CwdGuard guard;
  SetupPlasmaBasePath();

  DateTime dt = SampleEpoch();
  double r = 1.0 + 350.0 / RE;  // near the F2 peak
  std::vector<double> n = iri_2007(dt, r, 12.0, 0.2, 3.0);
  REQUIRE(n.size() == 4);
  for (double c : n) {
    REQUIRE(std::isfinite(c));
    REQUIRE(c >= 0.0);
  }
  REQUIRE(n[0] > 0.0);

  // Above 3000 km altitude the IRI helper returns zeros by design.
  double r_above = 1.0 + 4000.0 / RE;
  std::vector<double> n_above = iri_2007(dt, r_above, 12.0, 0.2, 3.0);
  REQUIRE(n_above.size() == 4);
  REQUIRE_THAT(n_above[0], WithinAbs(0.0, 1e-12));
}

// ----------------------------------- GCPM ----------------------------------

TEST_CASE("environment.plasma.gcpm.gcpm_v24_low_latitude_smoke") {
  CwdGuard guard;
  SetupPlasmaBasePath();
  set_iri_model(IRIModel::IRI_2007);  // the C++ ne_iri_ps_trough path needs an IRI model set

  DateTime dt = SampleEpoch();
  // A low magnetic latitude routes gcpm_v24 through ne_iri_ps_trough (the
  // plasmasphere/trough branch). The C++ path returns a finite, physical
  // electron density here. (The former NaN on this branch was not a translation
  // defect but a future-epoch solar-index data gap; see the future-epoch test
  // below and the year step-back fallback in iri_sm().)
  std::vector<double> out = gcpm_v24(dt, 1.0 + 320.0 / RE, 12.0, 0.2, 3.0);
  REQUIRE(out.size() == 8);
  INFO("gcpm_v24 (C++) low-latitude electron density: " << out[0]);
  REQUIRE(std::isfinite(out[0]));
  REQUIRE(out[0] > 0.0);
}

TEST_CASE("environment.plasma.gcpm.gcpm_v24_future_epoch") {
  CwdGuard guard;
  SetupPlasmaBasePath();
  set_iri_model(IRIModel::IRI_2007);

  // Epochs beyond the bundled solar-index projection (apf107.dat / ig_rz.dat end
  // ~2024, IRI projects only ~2 yr) used to make BOTH the C++ and Fortran GCPM
  // paths return NaN electron densities. The year step-back fallback (iri_sm and
  // the gcpm_v24_fortran wrapper) keeps them finite by falling back to the most
  // recent available solar activity.
  DateTime dt{2035, 172, 12, 0, 0.0};
  std::vector<double> out_cpp = gcpm_v24(dt, 1.0 + 320.0 / RE, 12.0, 0.2, 3.0);
  std::vector<double> out_fortran = gcpm_v24_fortran(dt, 1.0 + 320.0 / RE, 12.0, 0.2, 3.0);
  REQUIRE(std::isfinite(out_cpp[0]));
  REQUIRE(out_cpp[0] > 0.0);
  REQUIRE(std::isfinite(out_fortran[0]));
  REQUIRE(out_fortran[0] > 0.0);
}

TEST_CASE("environment.plasma.gcpm.gcpm_v24_fortran_density_profile") {
  CwdGuard guard;
  SetupPlasmaBasePath();

  DateTime dt = SampleEpoch();
  const double amlt = 12.0;
  const double alatr = 0.2;
  const double akp = 3.0;

  auto call = [&](double r) {
    std::vector<double> v = gcpm_v24_fortran(dt, r, amlt, alatr, akp);
    REQUIRE(v.size() == 4);
    for (double c : v) {
      REQUIRE(std::isfinite(c));
      REQUIRE(c >= 0.0);
    }
    return v;
  };

  // Ionospheric F-region (~320 km), mid-plasmasphere (2 Re) and outer
  // plasmasphere (4 Re). Total electron density must decrease monotonically
  // with geocentric distance.
  double ne_iono = call(1.0 + 320.0 / RE)[0];
  double ne_mid = call(2.0)[0];
  double ne_outer = call(4.0)[0];

  REQUIRE(ne_iono > 0.0);
  REQUIRE(ne_iono > ne_mid);
  REQUIRE(ne_mid > ne_outer);
  REQUIRE(ne_outer >= 0.0);
}

TEST_CASE("environment.plasma.gcpm.gcpm_v24_polar_cap") {
  CwdGuard guard;
  SetupPlasmaBasePath();

  DateTime dt = SampleEpoch();
  // High magnetic latitude drives a large L-shell -> polar-cap model branch.
  double alatr = 1.4;  // ~80 deg
  std::vector<double> out = gcpm_v24(dt, 1.0 + 350.0 / RE, 12.0, alatr, 3.0);
  REQUIRE(out.size() == 8);
  for (double c : out) REQUIRE(std::isfinite(c));
  REQUIRE(out[0] >= 0.0);
}

TEST_CASE("environment.plasma.gcpm.gcpm_v24_fortran") {
  CwdGuard guard;
  SetupPlasmaBasePath();

  DateTime dt = SampleEpoch();
  std::vector<double> out = gcpm_v24_fortran(dt, 1.0 + 320.0 / RE, 12.0, 0.2, 3.0);
  REQUIRE(out.size() == 4);
  for (double c : out) {
    REQUIRE(std::isfinite(c));
    REQUIRE(c >= 0.0);
  }
  // Total electron density in the daytime F-region is positive.
  REQUIRE(out[0] > 0.0);
}
