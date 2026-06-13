/**
 * @file ex_plasma_raytrace.cc
 * @author Keidai Iiyama
 * @brief Ray-trace GNSS links through the ionosphere/plasmasphere.
 * @version 0.1
 * @date 2025-02-17
 *
 * @copyright Copyright (c) 2025
 */

#include <omp.h>

#include <iostream>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;
using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {
  // input parameters
  if (argc < 10) {
    cerr << "Usage: " << argv[0]
         << " <sim_nums(int)> <cutoffRE(double)> <freq(int)> <stepsize(double)> "
            "<kp(double)> <use_moon(int)> <min_alt_low(double)> "
            "<min_alt_high(double)> <mainlobe_angle(double)>"
         << endl;
    return 1;  // Exit with error if no argument is provided
  }

  int sim_nums = stoi(argv[1]);           // Number of iterations to run
  double cutoff_RE = stod(argv[2]);       // Cutoff radius in Earth radii
  int freq_idx = stoi(argv[3]);           // GNSS frequency index: 1=L1, 2=L2, 5=L5
  double step_size = stod(argv[4]);       // Step size for ray tracing [km]
  double kp = stod(argv[5]);              // Kp index for the ionosphere model
  int use_moon = stoi(argv[6]);           // Use Moon in the simulation (1 for true, 0 for false)
  double min_alt_low = stod(argv[7]);     // Minimum altitude for visibility [km]
  double min_alt_high = stod(argv[8]);    // Maximum altitude for visibility [km]
  double mainlobe_angle = stod(argv[9]);  // Main lobe angle for visibility [degrees]  default is 20

  double freq;  // Frequency in Hz
  if (freq_idx == 1)
    freq = freq_L1;  // L1 frequency
  else if (freq_idx == 2)
    freq = freq_L2;  // L2 frequency
  else if (freq_idx == 5)
    freq = freq_L5;  // L5 frequency
  else {
    cerr << "Invalid frequency index. Use 1 for L1, 2 for L2, or 5 for L5." << endl;
    return 1;  // Exit with error if invalid frequency index is provided
  }

  // Simulation parameters ------------------------------------------------
  bool correction = true;              // Apply correction when ray tracing
  bool fine_correction = false;        // Apply fine correction when ray tracing
  double cutoff_r = cutoff_RE * RE;    // Cutoff radius for ray tracing [km]  (4 * RE)
  double margin_h = 100.0;             // Margin height for occultation [km]
  int num_Omega = 72;                  // Number of right ascension of ascending node points
  std::string integ_method = "Euler";  // Integration method ("RK4" or "Euler")
  std::string correction_method
      = "neldermead";       // Correction method ("grid" or "neldermead" "newton")
  bool debug_prop = false;  // Debug mode for propagation
  bool debug_corr = true;   // Debug mode for correction
  double corr_tol = 20.0;   // Correction tolerance [m] for ray tracing

  // setup the GNSS constellation ---------------------------------
  // The constellation file is resolved from the plasma data directory.
  std::vector<Satellite> gps_sats
      = setup_gnss_constellation("gps_2025_01_01.txt");  // Loading GPS constellation data

  // get the epoch
  double epoch_utc = gps_sats[0].epoch_utc_;
  double num_gps = gps_sats.size();

  // Setup a satellite in GEO
  int id_geo = 999;         // ID for GEO satellite
  double a_geo = 42164.0;   // Semi-major axis for GEO [km]
  double e_geo = 0.0001;    // Eccentricity for GEO (very small for circular orbit)
  double inc_geo = 0.0001;  // Inclination for GEO (close to 0 for equatorial orbit)

  // Moon
  double a_moon = 384400.0;         // Semi-major axis for Moon [km]
  double i_moon = 23.44 * DEG2RAD;  // Inclination for Moon [rad] (ecliptic plane)
  double e_moon = 0.0549;           // Eccentricity for Moon (slightly elliptical orbit)

  // Kp index for the ionosphere model
  DateTime datetime = mjd_to_datetime(tj2000_to_mjd(epoch_utc));  // Convert epoch to DateTime

  if (kp < 0) kp = get_kp_index(datetime);  // Get Kp index for the epoch

  // Create a RayTraceConfig object
  RayTraceConfig config;                     // Create a RayTraceConfig object
  config.freq_Hz = freq;                     // Set frequency in Hz
  config.step_size = step_size;              // Set step size for ray tracing [km]
  config.correction = correction;            // Apply correction to the ray tracing
  config.fine_correction = fine_correction;  // Apply fine correction to the ray
                                             // tracing
  config.cutoff_r = cutoff_r;                // Set cutoff radius for ray tracing [km]
  config.gradn_dx = 1.0;                     // Set gradient step size for refractive index [km]
  config.integ_method = integ_method;        // Set integration method ("
  config.kp = kp;                            // Kp index for the ionosphere model (default is -1)
  config.correction_method = correction_method;  // Set correction method
  config.use_fortran_gcpm = false;               // Use Fortran for ray tracing
  config.corr_tol = corr_tol;                    // Set correction tolerance [m] for ray tracing

  // set up the IRI model
  set_iri_model(IRIModel::IRI_2007);  // Set the IRI model to IRI-2007

  std::cout << "use Moon: " << (use_moon ? "Yes" : "No") << std::endl;
  std::cout << "Kp Index: " << kp << std::endl;                      // Print the Kp index
  std::cout << "frequency: " << freq * 1e-6 << " MHz" << std::endl;  // Print the frequency in MHz
  std::cout << " " << std::endl;

  std::vector<VecXd> store_vecs;  // Vector to store results
  int num_stored = 0;             // Number of stored vectors

  for (int omi = 0; omi < num_Omega; omi++) {
    double Omega
        = DEG2RAD * omi * (360 / num_Omega);  // Right ascension of ascending node for GEO [rad]

    // Create a receiver satellite in GEO or Moon
    double a_rx, e_rx, inc_rx;  // Semi-major axis, eccentricity, inclination
    if (use_moon) {
      a_rx = a_moon;    // Semi-major axis for Moon [km]
      e_rx = e_moon;    // Eccentricity for Moon
      inc_rx = i_moon;  // Inclination for Moon [rad]
    } else {
      a_rx = a_geo;      // Semi-major axis for GEO [km]
      e_rx = e_geo;      // Eccentricity for GEO
      inc_rx = inc_geo;  // Inclination for GEO [rad]
    }

    Satellite rx_sat(id_geo, Vec6d{a_rx, e_rx, inc_rx, Omega, 0.0, 0.0}, epoch_utc, GM_EARTH);

    for (int i = 0; i < num_gps; i++) {
      Vec3d pos_gps = gps_sats[i].get_pos();  // Get GPS satellite position
      Vec3d pos_rx = rx_sat.get_pos();

      // Compute occultation
      bool vis = compute_vis(pos_gps, pos_rx, RE + margin_h, mainlobe_angle);
      double min_alt = compute_min_altitude(pos_gps, pos_rx, RE);

      // trace the ray from the gps satellite to the receiver
      if ((vis) & (min_alt <= min_alt_high) & (min_alt >= min_alt_low)) {
        std::cout << num_stored << " | Omega:" << Omega * RAD2DEG << " ID: " << gps_sats[i].id_
                  << ", Minimum Altitude: " << min_alt << " km" << std::endl;

        Vec3d pos_tx = solve_lt(gps_sats[i], pos_rx,
                                epoch_utc);  // Solve light time for the transmission position

        VecXd store_vec(9);                // Vector to store the state derivative vector
        store_vec(0) = epoch_utc;          // Store epoch in seconds since J2000
        store_vec.segment<3>(1) = pos_tx;  // Store position of the transmitter
        store_vec.segment<3>(4) = pos_rx;  // Store position of the receiver
        store_vec(7) = gps_sats[i].id_;    // Store satellite
        store_vec(8) = min_alt;            // Store minimum altitude

        store_vecs.push_back(store_vec);  // Add the state derivative vector to
                                          // the vector of vectors
        num_stored++;                     // Increment the number of stored vectors
      }
    }  // End of loop over GPS satellites
  }  // End of loop over right ascension of ascending node points

  std::cout << "Total number of stored vectors: " << num_stored
            << std::endl;  // Print the total number of stored vectors
  std::cout << " " << std::endl;

  for (int i = 0; i < sim_nums; i++) {
    double comp_time = 0.0;  // Computation time for the ray tracing

    if (i >= num_stored) {
      std::cout << "No more stored vectors. Exiting." << std::endl;
      break;  // Exit if there are no more stored vectors
    }

    Vec3d pos_tx = store_vecs[i].segment<3>(1);  // Get position of the transmitter
    Vec3d pos_rx = store_vecs[i].segment<3>(4);  // Get position of the receiver

    std::cout << " =======================================================" << std::endl;
    std::cout << "  Raytrace ID: " << i << std::endl;
    std::cout << "  GPS ID: " << store_vecs[i](7) << std::endl;
    std::cout << "  Minimum Altitude: " << store_vecs[i](8) << " km" << std::endl;
    std::cout << " =========================================================" << std::endl;

    auto start_time = high_resolution_clock::now();  // Start timer
    PathProfile pp = trace_ray(epoch_utc, pos_tx, pos_rx, config, debug_prop, debug_corr);

    auto end_time = high_resolution_clock::now();  // End timer
    auto duration = duration_cast<seconds>(end_time - start_time);

    // print the total tec and delay
    std::cout << "  " << std::endl;
    std::cout << "[Raytrace Result] ID= " << i << std::endl;
    std::cout << "  TECU: " << pp.tecu << " TECU" << std::endl;
    std::cout << "  Total Delay     : " << pp.total_delay_m << " m" << std::endl;
    std::cout << "  Dist Total      : " << pp.sf << " m " << std::endl;
    std::cout << "  Dist Straight   : " << pp.dist_straight_km << " km" << std::endl;
    std::cout << "  Dist Bend       : " << pp.dist_bend_m << " m " << std::endl;
    std::cout << "  Total TEC Delay : " << pp.tec_delay_m << " m" << std::endl;
    std::cout << "  TEC Bend Delay  : " << pp.tec_delay_bend_m << " m" << std::endl;
    std::cout << "  Second Delay    : " << pp.second_delay_m << " m" << std::endl;
    std::cout << "  Final Pos Error : " << pp.corr_final_pos_err.norm() * 1000 << " m" << std::endl;
    std::cout << "  Final Time Error: " << pp.corr_final_time_err << " s" << std::endl;
    std::cout << "  Computation Time: " << duration.count() << " s" << std::endl;
    std::cout << "  " << std::endl;
  }

  return 0;
}
