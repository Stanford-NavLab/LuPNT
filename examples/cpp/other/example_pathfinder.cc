#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

#include <iostream>
#include <vector>

using namespace lupnt;

int main() {
  // Configuration
  bool show_orbit_plot = true;
  bool show_elevation_plot = true;

  // Epoch
  double epoch0 = (double)String2TAI("2025/10/02 00:00:00.000 UTC");

  // Orbital elements
  double sma = 5740;             // [km]  a, Semi-major axis
  double ecc = 0.58;             // [-]   e, Eccentricity
  double inc = deg2rad(54.856);  // [rad] i, Inclination
  double raan = deg2rad(0);  // [rad] W, Right ascension of the ascending node
  double aop = deg2rad(86.322);  // [rad] w, Argument of periapsis
  double ma = deg2rad(180);      // [rad] M, Mean anomaly
  Vec6 coe_sat_OP(sma, ecc, inc, raan, aop, ma);

  // Initial state
  auto rv_sat_OP = Classical2Cart(coe_sat_OP, GM_MOON);
  auto rv_sat_mi =
      ConvertFrame(epoch0, rv_sat_OP, Frame::MOON_OP, Frame::MOON_CI);

  // Time
  double T = 2 * M_PI * sqrt(pow(sma, 3) / GM_MOON);  // [s] Orbital period
  double dT = 0.5 * SECS_HOUR;  // [s] Data time step (2 hours)
  double dt = 5 * SECS_MINUTE;  // [s] Propagation time step (5 minutes)
  double tf = 2 * SECS_DAY;     // [s] Orbital period (14 days)
  int N_steps = (int)(tf / dT);

  // User positions
  double delta_lat = 5;                  // [deg]
  double delta_lon = 10;                 // [deg]
  double lat_min = -90, lat_max = 90;    // [deg]
  double lon_min = -180, lon_max = 180;  // [deg]

  int N_lats = (int)((lat_max - lat_min) / delta_lat + 1);
  int N_lons = (int)((lon_max - lon_min) / delta_lon + 1);
  auto lats = VecX::LinSpaced(N_lats, lat_min * RAD, lat_max * RAD);
  auto lons = VecX::LinSpaced(N_lons, lon_min * RAD, lon_max * RAD);
  MatX r_usr_pa(N_lats * N_lons, 3);  // [km] [x, y, z]
  for (int i = 0; i < N_lats; i++) {
    for (int j = 0; j < N_lons; j++) {
      int k = i * N_lons + j;
      Vec3 geo(lats(i), lons(j), 0);
      r_usr_pa.row(k) = LatLonAlt2Ecef(geo, R_MOON).transpose();
    }
  }
  std::cout << "Latitudes (deg)" << std::endl;
  std::cout << lats.transpose() * DEG << std::endl;
  std::cout << "Longitudes (deg)" << std::endl;
  std::cout << lons.transpose() * DEG << std::endl;

  // Dynamics
  NBodyDynamics dynamics;
  dynamics.SetPrimaryBody(Body::Moon());
  dynamics.AddBody(Body::Earth());
  dynamics.SetTimeStep(dt);

  // Main loop
  MatX rv_sat_mi_hist(N_steps, 6);  // [km,km/s] [x, y, z, vx, vy, vz]
  MatX rv_sat_pa_hist(N_steps, 6);  // [km,km/s] [x, y, z, vx, vy, vz]
  double epoch = epoch0;
  for (int i = 0; i < N_steps; i++) {
    // Propagate
    epoch += dT;
    dynamics.Propagate(rv_sat_mi, epoch - dT, epoch);

    auto rv_sat_pa =
        ConvertFrame(epoch, rv_sat_mi, Frame::MOON_CI, Frame::MOON_PA);
    rv_sat_mi_hist.row(i) = rv_sat_mi.transpose();
    rv_sat_pa_hist.row(i) = rv_sat_pa.transpose();
  }

  // Elevation and range
  MatX elevation(N_steps, N_lats * N_lons);  // [deg]
  MatX range(N_steps, N_lats * N_lons);      // [km]
  for (int i = 0; i < N_lats; i++) {
    for (int j = 0; j < N_lons; j++) {
      int k = i * N_lons + j;
      for (int t = 0; t < N_steps; t++) {
        Vec3 r_usr_pa_t = r_usr_pa.row(k).transpose();
        Vec3 r_sat_pa_t = rv_sat_pa_hist.row(t).head(3).transpose();
        auto [az, el, rng] = unpack(Cart2AzElRange(r_usr_pa_t, r_sat_pa_t));
        elevation(t, k) = rad2deg(el);
        range(t, k) = rng;
      }
    }
  }

  // Orbit plot
  double limit = 5 * R_MOON;
  VecX sizes = VecXd::Ones(N_lats * N_lons) * 1;

  if (show_orbit_plot) {
    matplot::figure();
    plot3(rv_sat_mi_hist.col(0), rv_sat_mi_hist.col(1), rv_sat_mi_hist.col(2));
    matplot::hold(matplot::on);
    matplot::xlabel("x [km]");
    matplot::ylabel("y [km]");
    matplot::zlabel("z [km]");
    matplot::xlim({-limit, limit});
    matplot::ylim({-limit, limit});
    matplot::zlim({-limit, limit});
    matplot::grid(matplot::on);
    scatter3(r_usr_pa.col(0), r_usr_pa.col(1), r_usr_pa.col(2), sizes);
    matplot::show();
  }

  // Elevation and range plot
  double min_elev = 5;                                         // [deg]
  VecX mask = VecX::Ones(N_steps) * NAN;                       // [-]
  VecX t = VecX::LinSpaced(N_steps, 0, tf);                    // [s]
  VecXd positions{{90, 0}, {60, -75}, {15, -90}, {-85, -30}};  // [deg]

  if (show_elevation_plot) {
    auto fig = matplot::figure();
    fig->position({0, 0, 1200, 800});
    fig->color("w");
    matplot::hold(matplot::on);
    for (int idx_pos = 0; idx_pos < positions.rows(); idx_pos++) {
      int idx_lat = (int)((90 - positions(idx_pos, 0)) / delta_lat);
      int idx_lon = (int)((positions(idx_pos, 1) + 180) / delta_lon);
      int idx_usr = (int)(idx_lat * N_lons + idx_lon);

      VecX elev = elevation.col(idx_usr);
      VecX rng = range.col(idx_usr);
      VecX elev_masked = (elev.array() >= min_elev).select(elev, mask);
      VecX rng_masked = (elev.array() >= min_elev).select(rng, mask);

      // Elevation
      matplot::subplot(positions.rows(), 2, 2 * idx_pos);
      matplot::hold(matplot::on);
      plot(t, elev_masked);
      matplot::plot(t, elev_masked);
      matplot::ylabel("Elevation [km]");
      matplot::xlabel("Time [s]");

      // Range
      matplot::subplot(positions.rows(), 2, 2 * idx_pos + 1);
      matplot::hold(matplot::on);
      plot(t, rng_masked);
      matplot::ylabel("Range [km]");
      matplot::xlabel("Time [s]");
    }
    matplot::show();
  }

  std::cout << "Done" << std::endl;
}