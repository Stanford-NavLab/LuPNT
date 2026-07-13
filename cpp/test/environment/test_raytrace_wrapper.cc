#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "lupnt/core/file.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/gcpm_interface.h"
#include "lupnt/environment/plasma/gcpm/iri_interface.h"
#include "lupnt/environment/plasma/tec/raytrace.h"
#include "lupnt/environment/plasma/tec/raytrace_wrapper.h"

using namespace pecsim;
using namespace Catch::Matchers;

// Covers the plasma TEC/ray-tracing entry points in raytrace.cc (get_iono_params,
// compute_B, refractive_index[_neB], gradient_n, compute_tec_straight, and the
// top-level trace_ray) plus the CSV loader in raytrace_wrapper.cc.
//
// Robustness: all electron-density evaluations use the Fortran GCPM path
// (use_fortran_gcpm = true), the physically-validated one. The pure-C++ GCPM has
// a known NaN branch (see test_plasma_gcpm.cc), so we assert finiteness /
// non-negativity rather than exact magnitudes. get_iono_params is the one entry
// point hard-wired to the C++ GCPM, so it is exercised at an in-range historical
// epoch where that path is finite.

namespace {
  struct CwdGuard {
    std::filesystem::path saved;
    CwdGuard() : saved(std::filesystem::current_path()) {}
    ~CwdGuard() { std::filesystem::current_path(saved); }
  };

  void SetupPlasma() {
    std::filesystem::path plasma_base = lupnt::GetDataPath() / "plasma";
    set_base_path(plasma_base.string());
    set_iri_model(IRIModel::IRI_2007);
  }

  // A daytime, mid-solar-cycle epoch inside the bundled apf107 / ig_rz range,
  // expressed as seconds since J2000.
  double SampleTj2000() {
    DateTime dt{2015, 172, 12, 0, 0.0};
    return mjd_to_tj2000(datetime_to_mjd(dt));
  }

  RayTraceConfig FortranConfig() {
    RayTraceConfig config;
    config.freq_Hz = freq_L1;
    config.kp = 3.0;
    config.use_fortran_gcpm = true;
    config.correction_method = "neldermead";  // default "grid" is rejected by trace_ray
    return config;
  }
}  // namespace

// ------------------------------ raytrace.cc --------------------------------

TEST_CASE("environment.plasma.raytrace.refractive_index_neB_pure") {
  const double freq = freq_L1;

  // With no plasma the phase refractive index is exactly 1.
  REQUIRE_THAT(refractive_index_neB(0.0, freq, 0.0, 0.0, false), WithinAbs(1.0, 1e-15));

  // A first-order (no magnetic) plasma lowers the phase index below 1, and more
  // electrons push it further down.
  double n1 = refractive_index_neB(1.0e11, freq, 0.0, 0.0, false);
  double n2 = refractive_index_neB(2.0e11, freq, 0.0, 0.0, false);
  REQUIRE(n1 < 1.0);
  REQUIRE(n2 < n1);
  REQUIRE(std::isfinite(n1));

  // Enabling the higher-order (magnetized) terms keeps the result finite and
  // still below unity for a positive electron density.
  double nq = refractive_index_neB(1.0e11, freq, 3.0e4, 0.5, true);
  REQUIRE(std::isfinite(nq));
  REQUIRE(nq < 1.0);
}

TEST_CASE("environment.plasma.raytrace.iono_params_and_B") {
  CwdGuard guard;
  SetupPlasma();
  double t = SampleTj2000();

  Vec3d params = get_iono_params(t, 3.0);
  for (int i = 0; i < 3; ++i) REQUIRE(std::isfinite(params[i]));
  REQUIRE(params[0] > 0.0);  // f107 solar flux
  REQUIRE(params[1] > 0.0);  // Rz12 sunspot number
  REQUIRE(params[2] > 0.0);  // hmF2 peak height [km]

  // Magnetic field just above the surface: finite and in a physical band
  // (~1e-5 .. 1e-4 in the units igrf14 returns; the plasma layer keeps them raw).
  RayTraceConfig config = FortranConfig();
  Vec3d pos(RE + 300.0, 0.0, 0.0);  // ECEF, km
  Vec3d B = compute_B(t, pos, config);
  for (int i = 0; i < 3; ++i) REQUIRE(std::isfinite(B[i]));
  REQUIRE(B.norm() > 0.0);
}

TEST_CASE("environment.plasma.raytrace.refractive_index_and_gradient") {
  CwdGuard guard;
  SetupPlasma();
  RayTraceConfig config = FortranConfig();
  double t = SampleTj2000();

  Vec3d pos(RE + 350.0, 0.0, 0.0);  // F-region, on the equatorial x-axis
  Vec3d shat(1.0, 0.0, 0.0);

  double n = refractive_index(t, pos, shat, config);
  REQUIRE(std::isfinite(n));
  // A daytime F-region plasma gives a phase index below unity but very close to it.
  REQUIRE(n < 1.0);
  REQUIRE(n > 0.9);

  Vec3d grad = gradient_n(t, pos, config);
  for (int i = 0; i < 3; ++i) REQUIRE(std::isfinite(grad[i]));
}

TEST_CASE("environment.plasma.raytrace.compute_tec_straight_nonnegative") {
  CwdGuard guard;
  SetupPlasma();
  RayTraceConfig config = FortranConfig();
  config.step_size = 100.0;  // coarse steps keep the integration quick
  double t = SampleTj2000();

  // A short ray climbing through the ionosphere on the equatorial x-axis.
  Vec3d tx(RE + 200.0, 0.0, 0.0);
  Vec3d rx(RE + 2000.0, 0.0, 0.0);
  double tec = compute_tec_straight(t, tx, rx, config);
  REQUIRE(std::isfinite(tec));
  REQUIRE(tec > 0.0);  // a daytime path accumulates positive electron content
}

TEST_CASE("environment.plasma.raytrace.trace_ray_no_correction") {
  CwdGuard guard;
  SetupPlasma();
  RayTraceConfig config = FortranConfig();
  config.correction = false;   // skip the Nelder-Mead bending correction
  config.straight_ray = true;  // straight-line propagation
  config.step_size = 100.0;
  double t_rx = SampleTj2000();

  // A downlink-style ray from high altitude down toward the ionosphere.
  Vec3d tx(RE + 20000.0, 0.0, 0.0);
  Vec3d rx(RE + 300.0, 0.0, 0.0);

  PathProfile pp = trace_ray(t_rx, tx, rx, config);
  REQUIRE(std::isfinite(pp.tecu));
  REQUIRE(pp.tecu >= 0.0);
  REQUIRE(std::isfinite(pp.tec_delay_m));
  REQUIRE(pp.tec_delay_m >= 0.0);
  REQUIRE(std::isfinite(pp.total_delay_m));
  REQUIRE(std::isfinite(pp.dist_straight_km));
  REQUIRE(pp.dist_straight_km > 0.0);

  // A ray whose reception epoch has no Kp configured is rejected downstream.
  RayTraceConfig bad = config;
  bad.kp = -1.0;
  REQUIRE_THROWS(trace_ray(t_rx, tx, rx, bad));
}

// -------------------------- raytrace_wrapper.cc ----------------------------

TEST_CASE("environment.plasma.raytrace.load_raytrace_data_csv") {
  namespace fs = std::filesystem;
  fs::path dir = fs::temp_directory_path() / "lupnt_raytrace_wrapper_test";
  fs::create_directories(dir);
  fs::path csv = dir / "rays.csv";

  {
    std::ofstream out(csv);
    out << "tidx,row_full,epoch_utc,min_alt,tx,rx,total_delay,tecu,tec_delay,"
           "second_delay,third_delay,dist_bend_m,tec_delay_bend_m,max_sep_line_m,"
           "final_pos_err_m\n";
    // A single well-formed row (15 comma-separated fields; vectors are
    // space-separated inside brackets so they survive the comma split).
    out << "2,5,1.5e8,100000,[7000000 0 0],[42000000 0 0],12.0,8.0,6.0,0.2,0.05,"
           "3.0,1.0,2.0,0.5\n";
    // A malformed row with too few fields is skipped.
    out << "1,2,3\n";
  }

  RaytraceData data = load_raytrace_data(csv.string());
  REQUIRE(data.num_rays() == 1);
  REQUIRE(data.tidxs.size() == 1);
  REQUIRE(data.tidxs[0] == 2);
  REQUIRE(data.row_fulls[0] == 5);
  REQUIRE_THAT(data.epoch_utcs[0], WithinRel(1.5e8, 1e-9));
  REQUIRE_THAT(data.min_alts[0], WithinAbs(100000.0, 1e-6));
  REQUIRE_THAT(data.tx_positions[0][0], WithinAbs(7000000.0, 1e-3));
  REQUIRE_THAT(data.rx_positions[0][0], WithinAbs(42000000.0, 1e-3));
  REQUIRE_THAT(data.total_delays[0], WithinAbs(12.0, 1e-9));
  REQUIRE_THAT(data.tecus[0], WithinAbs(8.0, 1e-9));
  REQUIRE_THAT(data.final_pos_err_m[0], WithinAbs(0.5, 1e-9));

  // A missing file yields an empty result rather than throwing.
  RaytraceData missing = load_raytrace_data((dir / "does_not_exist.csv").string());
  REQUIRE(missing.num_rays() == 0);

  fs::remove_all(dir);
}
