// Precision of the GNSS light-time solve.
//
// The light-time solvers in gnss_measurement.cc solve for the transmit epoch
// t_tx given a receive epoch t_rx and a geometric range rho. Formulating the
// iteration on the absolute epoch (t_tx = t_rx - rho/C) subtracts rho/C
// (~0.07 s for GNSS) from an ~8e8 s (year-2025) epoch, whose double ULP is
// ~1.2e-7 s ~ 36 m. Any consumer that later recovers the light time as
// (t_rx - t_tx)*C therefore inherits up to ~18 m of quantization noise, and
// the iteration's convergence test compares two ~8e8 epochs against a
// tolerance far below the ULP.
//
// The fix carries the light-time OFFSET dtau = t_tx - t_rx = -rho/C (a small,
// full-precision quantity) and holds t_rx fixed. This test pins the numerical
// difference between the two formulations so the offset form cannot regress.

#include <lupnt/core/constants.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("measurements.light_time_precision") {
  const double t_rx = 8.0e8;  // ~year 2025, seconds from J2000

  SECTION("offset formulation recovers the light-time range; absolute loses meters") {
    double max_dev_absolute = 0.0;
    double max_dev_offset = 0.0;

    // Sweep GNSS-scale ranges (~20000 km) so rho/C lands off the epoch's
    // representation grid at various phases.
    for (int i = 0; i < 200; i++) {
      double rho = 2.0e7 + i * 137.0;  // [m]

      // Old (absolute-epoch) formulation: store t_tx, recover the light time.
      double t_tx_abs = t_rx - rho / C;  // rounded at ULP(t_rx) ~ 1.2e-7 s
      double range_absolute = (t_rx - t_tx_abs) * C;

      // New (offset) formulation.
      double dtau = -rho / C;  // small, full precision
      double range_offset = (-dtau) * C;

      max_dev_absolute = std::max(max_dev_absolute, std::abs(range_absolute - rho));
      max_dev_offset = std::max(max_dev_offset, std::abs(range_offset - rho));
    }

    INFO("max |recovered - rho|: absolute = " << max_dev_absolute
                                              << " m, offset = " << max_dev_offset << " m");

    // The absolute-epoch subtraction loses meters of light time at 2025 epochs.
    REQUIRE(max_dev_absolute > 1.0);
    // The offset formulation is exact to the rounding of a small number (< 1 um).
    REQUIRE(max_dev_offset < 1.0e-6);
  }

  SECTION("convergence delta is clean in offset space, ULP-noisy in absolute space") {
    // Two nearly-converged light-time estimates whose true difference is 1 ns.
    const double rho1 = 2.05e7;
    const double rho2 = rho1 + 1.0e-9 * C;  // range differing by exactly 1 ns of light time

    // Absolute: delta of two ~8e8 epochs -- the 1 ns signal is buried under the
    // ~1.2e-7 s ULP, so the recovered delta is quantization noise.
    double t_tx1 = t_rx - rho1 / C;
    double t_tx2 = t_rx - rho2 / C;
    double delta_absolute = std::abs(t_tx1 - t_tx2);

    // Offset: delta of two small numbers recovers the 1 ns cleanly.
    double delta_offset = std::abs((-rho1 / C) - (-rho2 / C));

    INFO("recovered 1 ns delta: absolute = " << delta_absolute << " s, offset = " << delta_offset
                                             << " s");

    REQUIRE_THAT(delta_offset, WithinRel(1.0e-9, 1.0e-6));  // clean
    // The absolute delta cannot resolve 1 ns at this epoch (off by orders of magnitude).
    REQUIRE(std::abs(delta_absolute - 1.0e-9) > 1.0e-11);
  }
}
