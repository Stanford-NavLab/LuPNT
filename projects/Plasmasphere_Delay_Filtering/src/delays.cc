#include "src/delays.h"

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  Real ComputeShapiroDelay(double t_tai_rx, const Vec3& tx_pos, const Vec3& rx_pos, Frame frame) {
    // Constants (use your library constants if available)
    constexpr Real gamma = 1.0;  // GR
    const Real C3 = C * C * C;   // C = speed of light (m/s)
    const Real GM_sun = GM_SUN;  // [m^3/s^2] use your constant name

    // 1) Convert endpoints to GCRF if needed (positions in meters)
    Vec3 tx_gcrf
        = (frame == Frame::GCRF) ? tx_pos : ConvertFrame(t_tai_rx, tx_pos, frame, Frame::GCRF);
    Vec3 rx_gcrf
        = (frame == Frame::GCRF) ? rx_pos : ConvertFrame(t_tai_rx, rx_pos, frame, Frame::GCRF);

    // 2) First-cut light time -> approximate transmit epoch
    const Real R12 = (rx_gcrf - tx_gcrf).norm();
    const Real tau_geom = R12 / C;                // seconds
    const double t_tai_tx = t_tai_rx - tau_geom;  // simple 1-iteration improvement

    // 3) Get Sun position in the SAME frame as endpoints (here: geocentric Sun in GCRF)
    // You need something like: Sun position wrt Earth, expressed in GCRF at time t.
    // Replace this call with your own ephemeris interface.
    const Vec3 sun_gcrf_rx = GetBodyPos(t_tai_rx, BodyId::EARTH, BodyId::SUN, Frame::GCRF);
    const Vec3 sun_gcrf_tx = GetBodyPos(t_tai_tx, BodyId::EARTH, BodyId::SUN, Frame::GCRF);

    // 4) Distances from Sun to endpoints
    const Real r1 = (tx_gcrf - sun_gcrf_tx).norm();
    const Real r2 = (rx_gcrf - sun_gcrf_rx).norm();

    // 5) Guard against numerical issues: require r1 + r2 > R12
    const Real sum = r1 + r2;
    const Real den = sum - R12;
    const Real num = sum + R12;
    if (den <= Real(0) || num <= Real(0)) {
      return Real(0);  // geometry/inputs inconsistent; avoid log blow-up
    }

    // 6) Shapiro delay (seconds)
    const Real delay = ((1.0 + gamma) * GM_sun / C3) * log(num / den);
    return delay;
  }

  /**
   * @brief Compute Relativistic Delay due to Moon's gravity
   * @param rv_sat Satellite position and velocity [m] and [m/s]
   * @param frame Frame of the input position and velocity
   * @param dt Time since reference epoch [s] (since last epoch it was synced with LT)
   * @return Relativistic delay [s] with respect to LT (proper time - L)
   */
  Real ComputeRelativisticDelayLT(const Vec6& rv_sat, Frame frame, double dt) {
    // Compute Relativistic Delay from the Moon's gravitational field
    // rv_sat: satellite position and velocity in specified frame
    // dt: signal time of flight [s]
    Vec6 rv_sat_mci
        = (frame == Frame::MOON_CI) ? rv_sat : ConvertFrame(0.0, rv_sat, frame, Frame::MOON_CI);

    ClassicalOE coe = CartToClassical(Cart6(rv_sat_mci), GM_MOON);

    double C2 = C * C;

    double one_over_r0 = 0.0;

    bool use_LT = false;  // Set to true to use Lunar Time instead of LCT

    if (use_LT) {
      double u0 = 2.821e6;       // Moon gravitational potential at surface [m^2/s^2]
      double r0 = GM_MOON / u0;  // Reference radius [m]
      one_over_r0 = 1.0 / r0;
    }
    double E = MeanToEccAnomaly(coe.M(), coe.e());

    Real delay_1 = GM_MOON / C2 * (3 / 2 / coe.a() - one_over_r0) * dt;
    Real delay_2 = 2 / C2 * sqrt(GM_MOON * coe.a()) * coe.e() * sin(E);

    Real delay = delay_1 + delay_2;
    return delay;
  }

}  // namespace filtering_sim
