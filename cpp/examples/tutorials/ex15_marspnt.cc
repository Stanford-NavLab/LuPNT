// Example 15: Mars PNT -- a Walker navigation constellation around Mars
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex15_marspnt.ipynb.
//
// Builds a 9/3/1 Walker-delta navigation constellation at ~15,000 km altitude
// around Mars, propagates it for two Martian sidereal days with LuPNT's N-body
// dynamics (Mars 8x8 gravity field), and evaluates the resulting positioning
// service:
//
//   1. Acceleration budget on a satellite, decomposed by gravity term.
//   2. Walker constellation set-up in a Mars-equatorial frame, mapped to MARS_CI.
//   3. Propagation and conversion to the Mars body-fixed frame.
//   4. Geometric-dilution-of-precision (PDOP) and visibility over a lat/lon grid.
//   5. A single-point positioning fix at Jezero Crater (Perseverance site).
//
// The Python notebook renders the 3-D constellation, ground tracks and PDOP maps;
// this program prints the equivalent summary statistics.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "lupnt/lupnt.h"

using namespace lupnt;

namespace {

  // East-North-Up basis at a geodetic (lat, lon) on a sphere; rows are E, N, U.
  Eigen::Matrix3d EnuBasis(double lat, double lon) {
    const double cl = std::cos(lat), sl = std::sin(lat);
    const double co = std::cos(lon), so = std::sin(lon);
    Eigen::Matrix3d R;
    R.row(0) << -so, co, 0.0;
    R.row(1) << -sl * co, -sl * so, cl;
    R.row(2) << cl * co, cl * so, sl;
    return R;
  }

  // Classical Walker-delta elements [a,e,i,RAAN,argp,M] for T sats, P planes,
  // relative phasing F (columns match ClassicalToCart's COE layout).
  std::vector<Vec6> WalkerElements(int T, int P, int F, double a, double e, double inc) {
    const int spp = T / P;
    std::vector<Vec6> coes;
    for (int p = 0; p < P; ++p) {
      const double raan = p * 2.0 * M_PI / P;
      for (int s = 0; s < spp; ++s) {
        double M = s * 2.0 * M_PI / spp + p * 2.0 * M_PI * F / T;
        M = std::fmod(M, 2.0 * M_PI);
        coes.emplace_back(a, e, inc, raan, 0.0, M);
      }
    }
    return coes;
  }

}  // namespace

int main() {
  std::cout << std::fixed;

  // --- Mars constants from LuPNT's body model --------------------------------
  const Body mars_pm = Body::Mars();        // point-mass, for GM / R
  const double GM_MARS = mars_pm.GM.val();  // [m^3/s^2]
  const double R_MARS = mars_pm.R.val();    // [m]  mean volumetric radius

  // Mars sidereal rotation period, from the IAU prime-meridian rate LuPNT uses.
  const double W_DOT_MARS = 350.89198226;                     // [deg/day]
  const double SIDEREAL_DAY = 360.0 / W_DOT_MARS * SECS_DAY;  // [s]

  std::cout << std::setprecision(4);
  std::cout << "GM_MARS       = " << std::scientific << GM_MARS << " m^3/s^2\n" << std::fixed;
  std::cout << "R_MARS        = " << R_MARS / 1e3 << " km\n";
  std::cout << "Sidereal day  = " << SIDEREAL_DAY / 3600.0 << " h   (Mars: 24.6229 h)\n";

  // --- Force model: Mars 8x8 spherical-harmonic gravity ----------------------
  const int N_GRAV = 8;
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MARS_CI);  // Mars-centred inertial, ICRF-aligned
  dyn.SetUnits(SI_UNITS);
  dyn.AddBody(Body::Mars(N_GRAV, N_GRAV, "Mars50c.cof"));
  dyn.SetIntegrator(IntegratorType::RKF45);
  IntegratorParams params;
  params.reltol = 1e-11;
  params.abstol = 1e-3;  // [m]
  dyn.SetIntegratorParams(params);

  const Real t0 = GregorianToTime(2035, 1, 1, 0, 0, 0);

  // --- 1. Acceleration budget on a representative satellite ------------------
  Vec6 r_demo(R_MARS + 15000e3, 0.0, 0.0, 0.0, 2400.0, 2400.0);
  const std::map<std::string, Vec3> acc = dyn.ComputeAccelerations(t0, r_demo, true);
  std::vector<std::pair<std::string, double>> mags;
  for (const auto& [name, vec] : acc) mags.emplace_back(name, vec.norm().val());
  std::sort(mags.begin(), mags.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

  std::cout << "\nAcceleration budget at R_Mars + 15,000 km (Mars " << N_GRAV << "x" << N_GRAV
            << ", decomposed):\n";
  std::cout << std::scientific << std::setprecision(3);
  for (int i = 0; i < std::min<int>(8, mags.size()); ++i)
    std::cout << "  " << std::left << std::setw(16) << mags[i].first << std::right
              << "|a| = " << mags[i].second << " m/s^2\n";
  std::cout << std::fixed << std::setprecision(4);

  // --- 2. Mars-equatorial frame -> MARS_CI -----------------------------------
  // Rotation whose z-axis is the Mars pole (from the body-fixed frame), x-axis
  // the node of the Mars equator on the ICRF equator.
  const Vec3 pole_r
      = ConvertFrame(t0, Vec3(0.0, 0.0, 1.0), Frame::MARS_FIXED, Frame::MARS_CI, true);
  Eigen::Vector3d pole(pole_r.x().val(), pole_r.y().val(), pole_r.z().val());
  pole.normalize();
  Eigen::Vector3d xeq = Eigen::Vector3d::UnitZ().cross(pole).normalized();
  Eigen::Vector3d yeq = pole.cross(xeq);
  Eigen::Matrix3d R_eq2ci;
  R_eq2ci.col(0) = xeq;
  R_eq2ci.col(1) = yeq;
  R_eq2ci.col(2) = pole;
  std::cout << "Mars pole tilt from ICRF pole: " << std::acos(pole.z()) * DEG << " deg\n";

  // --- 3. Build the Walker constellation -------------------------------------
  const int N_SAT = 9, N_PLANES = 3, F_WALKER = 1;
  const int SATS_PER_PLANE = N_SAT / N_PLANES;
  const double ALT = 15000e3;
  const double A_SMA = R_MARS + ALT;
  const double INC = 45.0 * RAD;
  const double PERIOD = 2.0 * M_PI * std::sqrt(A_SMA * A_SMA * A_SMA / GM_MARS);

  const std::vector<Vec6> coes = WalkerElements(N_SAT, N_PLANES, F_WALKER, A_SMA, 0.0, INC);
  std::vector<Vec6> rv0(N_SAT);  // initial states in MARS_CI
  for (int k = 0; k < N_SAT; ++k) {
    const Vec6 rv_eq = ClassicalToCart(coes[k], GM_MARS);
    Eigen::Vector3d r_eq(rv_eq(0).val(), rv_eq(1).val(), rv_eq(2).val());
    Eigen::Vector3d v_eq(rv_eq(3).val(), rv_eq(4).val(), rv_eq(5).val());
    const Eigen::Vector3d r_ci = R_eq2ci * r_eq;
    const Eigen::Vector3d v_ci = R_eq2ci * v_eq;
    rv0[k] << r_ci.x(), r_ci.y(), r_ci.z(), v_ci.x(), v_ci.y(), v_ci.z();
  }
  std::cout << "\nWalker " << N_SAT << "/" << N_PLANES << "/" << F_WALKER
            << " constellation: " << N_PLANES << " planes x " << SATS_PER_PLANE
            << " sats, i = " << INC * DEG << " deg, alt " << ALT / 1e3 << " km\n";
  std::cout << "  orbital period " << PERIOD / 3600.0 << " h  vs Martian day "
            << SIDEREAL_DAY / 3600.0 << " h\n";

  // --- Propagate over two sols and rotate into the body-fixed frame ----------
  const double T_SIM = 2.0 * SIDEREAL_DAY;
  const double DT = 180.0;
  const VecX tspan = Arange(Real(0.0), Real(T_SIM), Real(DT));
  const VecX tfs = t0 + tspan.array();
  const int N_T = static_cast<int>(tfs.size());

  // pos_fixed[k] is an (N_T x 3) matrix of satellite positions in MARS_FIXED.
  std::vector<Eigen::MatrixXd> pos_fixed(N_SAT);
  for (int k = 0; k < N_SAT; ++k) {
    const MatX6 rv_ci = dyn.Propagate(rv0[k], tfs);
    const MatX3 r_fixed
        = ConvertFrame(tfs, MatX3(rv_ci.leftCols(3)), Frame::MARS_CI, Frame::MARS_FIXED, true);
    pos_fixed[k] = r_fixed.cast<double>();
  }
  std::cout << "Propagated " << N_SAT << " satellites over " << T_SIM / SIDEREAL_DAY << " sols ("
            << N_T << " epochs, " << T_SIM / PERIOD << " orbital periods)\n";

  // --- 4. Coverage and PDOP over a lat/lon grid ------------------------------
  const double MASK = 5.0 * RAD;
  std::vector<double> lat_g, lon_g;
  for (double la = -85.0; la <= 85.1; la += 5.0) lat_g.push_back(la);
  for (double lo = -180.0; lo < 180.0; lo += 5.0) lon_g.push_back(lo);

  std::vector<double> all_median_pdop;
  double vis_sum = 0.0;
  long vis_count = 0;
  int cells_below6 = 0, cells_with_fix = 0;
  for (double la_deg : lat_g) {
    for (double lo_deg : lon_g) {
      const double lat = la_deg * RAD, lon = lo_deg * RAD;
      const Eigen::Matrix3d Renu = EnuBasis(lat, lon);
      const Eigen::Vector3d r_user
          = R_MARS
            * Eigen::Vector3d(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon),
                              std::sin(lat));
      std::vector<double> pdop_t;
      for (int tt = 0; tt < N_T; ++tt) {
        std::vector<Eigen::Vector3d> los_enu;
        for (int k = 0; k < N_SAT; ++k) {
          const Eigen::Vector3d los = pos_fixed[k].row(tt).transpose() - r_user;
          const Eigen::Vector3d enu = Renu * los;
          const double rng = enu.norm();
          const double elev = std::asin(std::clamp(enu.z() / rng, -1.0, 1.0));
          if (elev >= MASK) los_enu.push_back(enu / rng);
        }
        vis_sum += los_enu.size();
        ++vis_count;
        if (static_cast<int>(los_enu.size()) < 4) continue;
        Eigen::MatrixXd H(los_enu.size(), 4);
        for (size_t s = 0; s < los_enu.size(); ++s) H.row(s) << los_enu[s].transpose(), 1.0;
        const Eigen::Matrix4d Q = (H.transpose() * H).inverse();
        const double pdop = std::sqrt(std::max(0.0, Q(0, 0) + Q(1, 1) + Q(2, 2)));
        pdop_t.push_back(pdop);
      }
      if (pdop_t.empty()) continue;
      std::sort(pdop_t.begin(), pdop_t.end());
      const double med = pdop_t[pdop_t.size() / 2];
      all_median_pdop.push_back(med);
      ++cells_with_fix;
      if (med < 6.0) ++cells_below6;
    }
  }
  std::sort(all_median_pdop.begin(), all_median_pdop.end());
  const double global_median_pdop
      = all_median_pdop.empty() ? std::nan("") : all_median_pdop[all_median_pdop.size() / 2];
  std::cout << "\nService volume (" << lat_g.size() << "x" << lon_g.size()
            << " lat/lon grid, 5 deg mask):\n";
  std::cout << "  global median PDOP               : " << global_median_pdop << "\n";
  std::cout << "  grid fraction with median PDOP<6 : "
            << 100.0 * cells_below6 / std::max(1, cells_with_fix) << " %\n";
  std::cout << "  mean satellites in view          : " << vis_sum / std::max(1L, vis_count) << "\n";

  // --- 5. Single-point positioning at Jezero Crater --------------------------
  const double JEZERO_LAT = 18.44 * RAD, JEZERO_LON = 77.45 * RAD;
  const Eigen::Matrix3d Renu = EnuBasis(JEZERO_LAT, JEZERO_LON);
  const Eigen::Vector3d r_user
      = R_MARS
        * Eigen::Vector3d(std::cos(JEZERO_LAT) * std::cos(JEZERO_LON),
                          std::cos(JEZERO_LAT) * std::sin(JEZERO_LON), std::sin(JEZERO_LAT));
  const double SIGMA = 8.0;  // [m] pseudorange noise
  RandomEngine::SetSeed(42);

  std::vector<double> err3d, err_h, err_v, pdop_fix;
  for (int tt = 0; tt < N_T; ++tt) {
    std::vector<Eigen::Vector3d> sats;
    for (int k = 0; k < N_SAT; ++k) {
      const Eigen::Vector3d s = pos_fixed[k].row(tt).transpose();
      const Eigen::Vector3d enu = Renu * (s - r_user);
      const double elev = std::asin(std::clamp(enu.z() / enu.norm(), -1.0, 1.0));
      if (elev >= MASK) sats.push_back(s);
    }
    const int m = static_cast<int>(sats.size());
    if (m < 4) continue;

    // Synthetic pseudoranges (true range + white noise).
    Eigen::VectorXd rho(m);
    for (int s = 0; s < m; ++s)
      rho(s) = (sats[s] - r_user).norm() + SIGMA * SampleNormal(0.0, 1.0).val();

    // Gauss-Newton single-point fix for [x, y, z, c*dt].
    Eigen::Vector4d x = Eigen::Vector4d::Zero();
    Eigen::MatrixXd H(m, 4);
    for (int it = 0; it < 15; ++it) {
      Eigen::VectorXd res(m);
      for (int s = 0; s < m; ++s) {
        const Eigen::Vector3d d = x.head(3) - sats[s];
        const double r = d.norm();
        H.row(s) << (d / r).transpose(), 1.0;
        res(s) = rho(s) - (r + x(3));
      }
      const Eigen::Vector4d dx = (H.transpose() * H).ldlt().solve(H.transpose() * res);
      x += dx;
      if (dx.head(3).norm() < 1e-3) break;
    }
    const Eigen::Vector3d derr = Renu * (x.head(3) - r_user);
    err3d.push_back(derr.norm());
    err_h.push_back(std::hypot(derr.x(), derr.y()));
    err_v.push_back(std::abs(derr.z()));
    const Eigen::Matrix4d Q = (H.transpose() * H).inverse();
    pdop_fix.push_back(std::sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2)));
  }

  auto median = [](std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v.empty() ? std::nan("") : v[v.size() / 2];
  };
  auto pct = [](std::vector<double> v, double p) {
    std::sort(v.begin(), v.end());
    return v.empty() ? std::nan("")
                     : v[std::min<size_t>(v.size() - 1, size_t(p / 100.0 * v.size()))];
  };
  std::cout << "\nPositioning at Jezero Crater (18.44 N, 77.45 E), " << SIGMA
            << " m pseudorange noise:\n";
  std::cout << "  epochs with >=4 satellites : " << err3d.size() << " / " << N_T << " ("
            << 100.0 * err3d.size() / N_T << " %)\n";
  std::cout << "  3-D error   median " << median(err3d) << " m   95% " << pct(err3d, 95) << " m\n";
  std::cout << "    horizontal median " << median(err_h) << " m\n";
  std::cout << "    vertical   median " << median(err_v) << " m\n";
  std::cout << "  PDOP at fix epochs median " << median(pdop_fix) << "\n";

  return 0;
}
