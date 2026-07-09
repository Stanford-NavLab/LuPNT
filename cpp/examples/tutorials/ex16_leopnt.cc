// Example 16: LEO PNT -- a low-Earth-orbit navigation constellation
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex16_leopnt.ipynb.
//
// The Earth companion to the Mars tutorial (Example 15). Builds a dense low-Earth
// -orbit Position/Navigation/Timing constellation -- the kind proposed to augment
// GPS with strong, close-range signals -- and evaluates the positioning service it
// delivers to a user in San Francisco. LEO adds a force higher orbits can ignore:
// atmospheric drag, which LuPNT models with a Harris-Priester density model.
//
//   1. Force model: J2 (EGM96) gravity + Harris-Priester drag; acceleration budget.
//   2. Why drag matters: unmodeled-drag ephemeris drift over five days.
//   3. A 110/10/1 Walker constellation at 600 km, 55 deg, propagated with J2 + drag.
//   4. Coverage and PDOP over a lat/lon grid.
//   5. Single-point positioning at San Francisco.
//   6. On-board orbit- and time-determination (ODTS) with GPS via an EKF.
//
// The Python notebook renders the 3-D constellation, ground tracks, PDOP maps and
// filter-convergence plots; this program prints the equivalent summary statistics.

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

// True if the LEO->GPS line of sight clears the Earth (occulting radius r_occ).
bool GpsVisible(const Eigen::Vector3d& r_leo, const Eigen::Vector3d& r_gps, double r_occ) {
  const Eigen::Vector3d d = r_gps - r_leo;
  const double L = d.norm();
  const Eigen::Vector3d u = d / L;
  const double s = -r_leo.dot(u);  // foot of perpendicular from Earth centre
  if (s < 0.0 || s > L) return true;
  return (r_leo + s * u).norm() > r_occ;
}

double Median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v.empty() ? std::nan("") : v[v.size() / 2];
}
double Percentile(std::vector<double> v, double p) {
  std::sort(v.begin(), v.end());
  return v.empty() ? std::nan("")
                   : v[std::min<size_t>(v.size() - 1, size_t(p / 100.0 * v.size()))];
}

}  // namespace

int main() {
  std::cout << std::fixed;

  // --- 0. Earth constants from LuPNT's body model ----------------------------
  const Body earth_pm = Body::Earth();      // point-mass, for GM / R
  const double GM_EARTH = earth_pm.GM.val();
  const double R_EARTH = earth_pm.R.val();
  std::cout << std::setprecision(4);
  std::cout << "GM_EARTH      = " << std::scientific << GM_EARTH << " m^3/s^2\n" << std::fixed;
  std::cout << "R_EARTH       = " << R_EARTH / 1e3 << " km\n";

  // --- 1. Force model: J2 (EGM96 4x4) gravity + Harris-Priester drag ---------
  const int N_GRAV = 4;
  const double CD = 2.2, AREA = 1.0, MASS = 100.0;
  const double BSTAR = CD * AREA / MASS;  // ballistic coefficient B = CD*A/m [m^2/kg]
  const Body earth = Body::Earth(N_GRAV, N_GRAV, "EGM96.cof");
  const Frame GCRF = earth.inertial_frame;
  const Frame ITRF = earth.fixed_frame;

  auto make_dynamics = [&](bool drag, IntegratorType integ, double step) {
    NBodyDynamics dyn;
    dyn.SetFrame(GCRF);
    dyn.SetUnits(SI_UNITS);
    dyn.AddBody(earth);
    dyn.SetIntegrator(integ);
    dyn.SetTimeStep(step);
    if (drag) dyn.SetDragCoeff(BSTAR);
    return dyn;
  };

  NBodyDynamics dyn = make_dynamics(true, IntegratorType::RK4, 30.0);
  const Real EPOCH = GregorianToTime(2035, 1, 1, 0, 0, 0);

  std::cout << "\nBallistic coefficient B = CD*A/m = " << std::setprecision(3) << BSTAR
            << " m^2/kg\n";
  std::cout << std::setprecision(3) << std::scientific;
  std::cout << std::setw(10) << "altitude" << std::setw(14) << "|a_grav|" << std::setw(14)
            << "|a_J2|" << std::setw(14) << "|a_drag|" << "\n";
  for (double alt : {400e3, 600e3, 800e3}) {
    Vec6 coe(R_EARTH + alt, 0.001, 55.0 * RAD, 0.0, 0.0, 0.0);
    Vec6 rv = ClassicalToCart(coe, GM_EARTH);
    const std::map<std::string, Vec3> acc = dyn.ComputeAccelerations(EPOCH, rv, true);
    auto mag = [&](const std::string& k) {
      auto it = acc.find(k);
      return it == acc.end() ? 0.0 : it->second.norm().val();
    };
    std::cout << std::fixed << std::setprecision(0) << std::setw(7) << alt / 1e3 << " km"
              << std::scientific << std::setprecision(3) << std::setw(14) << mag("EARTH_gravity")
              << std::setw(14) << mag("EARTH_J2") << std::setw(14) << mag("drag") << "\n";
  }
  std::cout << std::fixed << std::setprecision(4);

  // --- 2. Why drag matters: unmodeled-drag ephemeris drift -------------------
  const double ALT_DEMO = 600e3;
  Vec6 coe_demo(R_EARTH + ALT_DEMO, 0.0005, 55.0 * RAD, 0.0, 0.0, 0.0);
  Vec6 rv_demo = ClassicalToCart(coe_demo, GM_EARTH);
  const VecX t_demo = EPOCH + Arange(Real(0.0), Real(5.0 * 86400.0), Real(900.0)).array();
  NBodyDynamics dyn_drag = make_dynamics(true, IntegratorType::RK8, 60.0);
  NBodyDynamics dyn_nodrag = make_dynamics(false, IntegratorType::RK8, 60.0);
  const MatX6 xd = dyn_drag.Propagate(rv_demo, t_demo);
  const MatX6 xn = dyn_nodrag.Propagate(rv_demo, t_demo);
  double sep_last = 0.0;
  for (int i = 0; i < t_demo.size(); ++i)
    sep_last = (xd.row(i).head(3) - xn.row(i).head(3)).norm().val();
  std::cout << "\nNeglecting drag mispredicts the 600 km satellite by " << sep_last / 1e3
            << " km after 5 days.\n";

  // --- 3. A 110/10/1 Walker constellation at 600 km --------------------------
  const int N_SAT = 110, N_PLANES = 10, F_WALKER = 1;
  const int SATS_PER_PLANE = N_SAT / N_PLANES;
  const double ALT = 600e3;
  const double A_SMA = R_EARTH + ALT;
  const double INC = 55.0 * RAD;
  const double PERIOD = 2.0 * M_PI * std::sqrt(A_SMA * A_SMA * A_SMA / GM_EARTH);
  const std::vector<Vec6> coes = WalkerElements(N_SAT, N_PLANES, F_WALKER, A_SMA, 0.001, INC);
  std::vector<Vec6> rv0(N_SAT);
  for (int k = 0; k < N_SAT; ++k) rv0[k] = ClassicalToCart(coes[k], GM_EARTH);
  std::cout << "\nWalker " << N_SAT << "/" << N_PLANES << "/" << F_WALKER << ": " << N_PLANES
            << " planes x " << SATS_PER_PLANE << " sats, i = " << INC * DEG << " deg, alt "
            << ALT / 1e3 << " km\n";
  std::cout << "  orbital period " << PERIOD / 60.0 << " min (" << 86400.0 / PERIOD
            << " revs/day)\n";

  // Propagate the fleet over two orbital periods and rotate GCRF -> ITRF.
  const double T_SIM = 2.0 * PERIOD;
  const double DT = 60.0;
  const VecX TS = EPOCH + Arange(Real(0.0), Real(T_SIM), Real(DT)).array();
  const int N_T = static_cast<int>(TS.size());
  std::vector<Eigen::MatrixXd> pos_itrf(N_SAT);
  for (int k = 0; k < N_SAT; ++k) {
    const MatX6 rv_gcrf = dyn.Propagate(rv0[k], TS);
    const MatX3 r_itrf = ConvertFrame(TS, MatX3(rv_gcrf.leftCols(3)), GCRF, ITRF, true);
    pos_itrf[k] = r_itrf.cast<double>();
  }
  std::cout << "Propagated " << N_SAT << " satellites over " << T_SIM / 60.0 << " min (" << N_T
            << " epochs) with J2 + drag\n";

  // --- 4. Coverage and PDOP over a lat/lon grid ------------------------------
  const double MASK = 5.0 * RAD;
  std::vector<double> lat_g, lon_g;
  for (double la = -82.0; la <= 82.1; la += 4.0) lat_g.push_back(la);
  for (double lo = -180.0; lo < 180.0; lo += 4.0) lon_g.push_back(lo);
  std::vector<double> band_median_pdop;
  double band_vis_sum = 0.0;
  long band_vis_count = 0;
  for (double la_deg : lat_g) {
    for (double lo_deg : lon_g) {
      const double lat = la_deg * RAD, lon = lo_deg * RAD;
      const Eigen::Matrix3d Renu = EnuBasis(lat, lon);
      const Eigen::Vector3d r_user =
          R_EARTH * Eigen::Vector3d(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon),
                                    std::sin(lat));
      const bool in_band = std::abs(la_deg) <= 60.0;
      std::vector<double> pdop_t;
      for (int tt = 0; tt < N_T; tt += 2) {  // subsample time to keep the map fast
        std::vector<Eigen::Vector3d> los_enu;
        for (int k = 0; k < N_SAT; ++k) {
          const Eigen::Vector3d enu = Renu * (pos_itrf[k].row(tt).transpose() - r_user);
          const double rng = enu.norm();
          if (std::asin(std::clamp(enu.z() / rng, -1.0, 1.0)) >= MASK) los_enu.push_back(enu / rng);
        }
        if (in_band) {
          band_vis_sum += los_enu.size();
          ++band_vis_count;
        }
        if (static_cast<int>(los_enu.size()) < 4) continue;
        Eigen::MatrixXd H(los_enu.size(), 4);
        for (size_t s = 0; s < los_enu.size(); ++s) H.row(s) << los_enu[s].transpose(), 1.0;
        const Eigen::Matrix4d Q = (H.transpose() * H).inverse();
        pdop_t.push_back(std::sqrt(std::max(0.0, Q(0, 0) + Q(1, 1) + Q(2, 2))));
      }
      if (in_band && !pdop_t.empty()) {
        std::sort(pdop_t.begin(), pdop_t.end());
        band_median_pdop.push_back(pdop_t[pdop_t.size() / 2]);
      }
    }
  }
  std::cout << "\nService band (|lat| <= 60 deg, 5 deg mask):\n";
  std::cout << "  median PDOP              : " << Median(band_median_pdop) << "\n";
  std::cout << "  mean satellites in view  : " << band_vis_sum / std::max(1L, band_vis_count)
            << "\n";

  // --- 5. Single-point positioning at San Francisco --------------------------
  const double SF_LAT = 37.7749 * RAD, SF_LON = -122.4194 * RAD;
  const Eigen::Matrix3d Renu = EnuBasis(SF_LAT, SF_LON);
  const Eigen::Vector3d r_user =
      R_EARTH * Eigen::Vector3d(std::cos(SF_LAT) * std::cos(SF_LON),
                                std::cos(SF_LAT) * std::sin(SF_LON), std::sin(SF_LAT));
  const double SIGMA = 1.0;  // [m] ranging noise
  RandomEngine::SetSeed(0);
  std::vector<double> err3d, err_h, err_v, pdop_fix;
  int n_vis_min = N_SAT, n_vis_max = 0;
  double n_vis_sum = 0.0;
  for (int tt = 0; tt < N_T; ++tt) {
    std::vector<Eigen::Vector3d> sats;
    for (int k = 0; k < N_SAT; ++k) {
      const Eigen::Vector3d s = pos_itrf[k].row(tt).transpose();
      const Eigen::Vector3d enu = Renu * (s - r_user);
      if (std::asin(std::clamp(enu.z() / enu.norm(), -1.0, 1.0)) >= MASK) sats.push_back(s);
    }
    const int m = static_cast<int>(sats.size());
    n_vis_min = std::min(n_vis_min, m);
    n_vis_max = std::max(n_vis_max, m);
    n_vis_sum += m;
    if (m < 4) continue;
    Eigen::VectorXd rho(m);
    Eigen::Vector3d sat_mean = Eigen::Vector3d::Zero();
    for (int s = 0; s < m; ++s) {
      rho(s) = (sats[s] - r_user).norm() + SIGMA * SampleNormal(0.0, 1.0).val();
      sat_mean += sats[s];
    }
    sat_mean /= m;
    // Gauss-Newton single-point fix for [x, y, z, c*dt]; start on the surface
    // beneath the mean satellite direction (the notebook's initial guess).
    Eigen::Vector4d x = Eigen::Vector4d::Zero();
    x.head(3) = R_EARTH * sat_mean.normalized();
    Eigen::MatrixXd H(m, 4);
    bool converged = false;
    for (int it = 0; it < 20; ++it) {
      Eigen::VectorXd res(m);
      for (int s = 0; s < m; ++s) {
        const Eigen::Vector3d d = x.head(3) - sats[s];
        const double r = d.norm();
        H.row(s) << (d / r).transpose(), 1.0;
        res(s) = rho(s) - (r + x(3));
      }
      const Eigen::Vector4d dx = (H.transpose() * H).ldlt().solve(H.transpose() * res);
      x += dx;
      if (dx.head(3).norm() < 1e-3) {
        converged = true;
        break;
      }
    }
    if (!converged || std::abs(x.head(3).norm() - R_EARTH) > 5e5) continue;
    const Eigen::Vector3d derr = Renu * (x.head(3) - r_user);
    err3d.push_back(derr.norm());
    err_h.push_back(std::hypot(derr.x(), derr.y()));
    err_v.push_back(std::abs(derr.z()));
    const Eigen::Matrix4d Q = (H.transpose() * H).inverse();
    pdop_fix.push_back(std::sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2)));
  }
  std::cout << "\nPositioning at San Francisco (37.77 N, 122.42 W), " << SIGMA
            << " m ranging noise:\n";
  std::cout << "  satellites in view : min " << n_vis_min << ", mean " << n_vis_sum / N_T
            << ", max " << n_vis_max << "\n";
  std::cout << "  epochs with a fix  : " << err3d.size() << " / " << N_T << " ("
            << 100.0 * err3d.size() / N_T << " %)\n";
  std::cout << "  3-D error median " << Median(err3d) << " m   95% " << Percentile(err3d, 95)
            << " m\n";
  std::cout << "    horizontal median " << Median(err_h) << " m\n";
  std::cout << "    vertical   median " << Median(err_v) << " m\n";
  std::cout << "  PDOP at fix epochs median " << Median(pdop_fix) << "\n";

  // --- 6. On-board orbit- and time-determination (ODTS) with GPS -------------
  // A nominal GPS Walker 24/6/2 at ~20,200 km serves as the truth reference; the
  // LEO nav-satellite carries a GPS receiver and runs an EKF over [r, v, b, d].
  const std::vector<Vec6> gps_coes =
      WalkerElements(24, 6, 2, R_EARTH + 20200e3, 0.0, 55.0 * RAD);
  CartesianTwoBodyDynamics gps_dyn(GM_EARTH);
  const double DT_ODTS = 30.0;
  const VecX TS_ODTS = EPOCH + Arange(Real(0.0), Real(2.0 * PERIOD), Real(DT_ODTS)).array();
  const int N_ODTS = static_cast<int>(TS_ODTS.size());

  // Truth: LEO satellite 0 WITH drag, and the GPS satellites (two-body).
  const MatX6 leo_truth = dyn.Propagate(rv0[0], TS_ODTS);
  std::vector<Eigen::MatrixXd> gps_truth(gps_coes.size());
  for (size_t j = 0; j < gps_coes.size(); ++j) {
    const MatX6 g = gps_dyn.Propagate(ClassicalToCart(gps_coes[j], GM_EARTH), TS_ODTS);
    gps_truth[j] = g.leftCols(3).cast<double>();
  }
  const double R_OCC = R_EARTH + 50e3;

  NBodyDynamics dyn_filter = make_dynamics(false, IntegratorType::RK4, 30.0);  // gravity only
  const double SIGMA_RHO = 1.0, Q_ACC = 1e-5, Q_CLK = 5e-4;
  RandomEngine::SetSeed(1);

  // Truth receiver clock: bias b [m] and drift d [m/s] (random walk).
  std::vector<double> b_true(N_ODTS, 0.0), d_true(N_ODTS, 0.0);
  b_true[0] = 30.0;
  d_true[0] = 0.20;
  for (int k = 1; k < N_ODTS; ++k) {
    d_true[k] = d_true[k - 1] + SampleNormal(0.0, 1.0).val() * 1e-4 * std::sqrt(DT_ODTS);
    b_true[k] = b_true[k - 1] + d_true[k - 1] * DT_ODTS;
  }

  // Initial estimate: truth + gross errors.
  Eigen::VectorXd x(8);
  for (int i = 0; i < 6; ++i) x(i) = rv0[0](i).val();
  x(6) = b_true[0] + 80.0;
  x(7) = d_true[0] + 0.4;
  for (int i = 0; i < 3; ++i) x(i) += SampleNormal(0.0, 1.0).val() * 2000.0;
  for (int i = 3; i < 6; ++i) x(i) += SampleNormal(0.0, 1.0).val() * 2.0;
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(8, 8);
  P.diagonal() << 3e3 * 3e3, 3e3 * 3e3, 3e3 * 3e3, 3.0, 3.0, 3.0, 300.0 * 300.0, 1.0;

  std::vector<double> pos_err, pos_3sig, clk_err;
  std::vector<int> gps_in_view;
  for (int k = 0; k < N_ODTS; ++k) {
    if (k > 0) {
      // Predict: native STM for the orbit (autodiff), analytic model for the clock.
      Vec6 x_orbit;
      for (int i = 0; i < 6; ++i) x_orbit(i) = x(i);
      MatXd stm(6, 6);
      TerminationInfo info;
      const State xf = dyn_filter.PropagateExStm(x_orbit, TS_ODTS[k - 1], TS_ODTS[k], &stm, &info);
      const VecXd xf_v = xf.cast<double>();
      for (int i = 0; i < 6; ++i) x(i) = xf_v(i);
      x(6) += x(7) * DT_ODTS;
      Eigen::MatrixXd Phi = Eigen::MatrixXd::Identity(8, 8);
      Phi.block(0, 0, 6, 6) = stm;
      Phi(6, 7) = DT_ODTS;
      Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(8, 8);
      const double dt = DT_ODTS;
      Q.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity() * Q_ACC * Q_ACC * dt * dt * dt / 3.0;
      Q.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity() * Q_ACC * Q_ACC * dt * dt / 2.0;
      Q.block(3, 0, 3, 3) = Q.block(0, 3, 3, 3).transpose();
      Q.block(3, 3, 3, 3) = Eigen::Matrix3d::Identity() * Q_ACC * Q_ACC * dt;
      Q(6, 6) = Q_CLK * Q_CLK * dt * dt * dt / 3.0;
      Q(7, 7) = Q_CLK * Q_CLK * dt;
      P = Phi * P * Phi.transpose() + Q;
    }

    // A-priori errors (before this epoch's measurement update).
    const Eigen::Vector3d r_true = leo_truth.row(k).head(3).cast<double>().transpose();
    pos_err.push_back((x.head(3) - r_true).norm());
    pos_3sig.push_back(3.0 * std::sqrt(P(0, 0) + P(1, 1) + P(2, 2)));
    clk_err.push_back(x(6) - b_true[k]);

    // Update: one pseudorange per visible GPS satellite.
    std::vector<int> vis;
    for (size_t j = 0; j < gps_coes.size(); ++j)
      if (GpsVisible(r_true, gps_truth[j].row(k).transpose(), R_OCC)) vis.push_back(j);
    gps_in_view.push_back(static_cast<int>(vis.size()));
    const int mv = static_cast<int>(vis.size());
    if (mv == 0) continue;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(mv, 8);
    Eigen::VectorXd y(mv);
    for (int i = 0; i < mv; ++i) {
      const Eigen::Vector3d r_gps = gps_truth[vis[i]].row(k).transpose();
      const double rho = (r_gps - r_true).norm() + b_true[k] + SampleNormal(0.0, 1.0).val() * SIGMA_RHO;
      const double rho_pred = (r_gps - x.head(3)).norm();
      H.block(i, 0, 1, 3) = (-(r_gps - x.head(3)) / rho_pred).transpose();
      H(i, 6) = 1.0;
      y(i) = rho - (rho_pred + x(6));
    }
    const Eigen::MatrixXd S =
        H * P * H.transpose() + Eigen::MatrixXd::Identity(mv, mv) * SIGMA_RHO * SIGMA_RHO;
    const Eigen::MatrixXd K = P * H.transpose() * S.inverse();
    x += K * y;
    P = (Eigen::MatrixXd::Identity(8, 8) - K * H) * P;
  }

  const int half = N_ODTS / 2;
  auto rms = [&](const std::vector<double>& v) {
    double s = 0.0;
    for (int i = half; i < N_ODTS; ++i) s += v[i] * v[i];
    return std::sqrt(s / (N_ODTS - half));
  };
  int gmin = 99, gmax = 0;
  double gsum = 0.0;
  for (int g : gps_in_view) {
    gmin = std::min(gmin, g);
    gmax = std::max(gmax, g);
    gsum += g;
  }
  int consistent = 0;
  for (int i = half; i < N_ODTS; ++i)
    if (pos_err[i] < pos_3sig[i]) ++consistent;
  std::cout << "\nOn-board ODTS with GPS (EKF over [r, v, clock bias, drift]):\n";
  std::cout << "  GPS in view from LEO : min " << gmin << ", mean " << gsum / N_ODTS << ", max "
            << gmax << "\n";
  std::cout << "  position: initial error " << pos_err[0] / 1e3 << " km  ->  steady-state RMS "
            << rms(pos_err) << " m\n";
  std::cout << "  clock bias: steady-state estimation-error RMS " << rms(clk_err) << " m ("
            << rms(clk_err) / 299792458.0 * 1e9 << " ns)\n";
  std::cout << "  filter consistency: " << 100.0 * consistent / (N_ODTS - half)
            << " % of steady-state epochs have |error| < 3-sigma\n";

  return 0;
}
