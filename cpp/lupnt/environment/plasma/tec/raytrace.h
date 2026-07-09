/**
 * @file tec.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-17
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <functional>
#include <limits>
#include <vector>

#include "lupnt/environment/plasma/core/definitions.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"
#include "lupnt/environment/plasma/gcpm/conversions.h"
#include "lupnt/environment/plasma/gcpm/gcpm_interface.h"

namespace pecsim {

  // Constants
  const double freq_L1 = 1575.42e6;  // L1 frequency
  const double freq_L2 = 1227.60e6;  // L2 frequency
  const double freq_L5 = 1176.45e6;  // L5 frequency

  struct RayTraceConfig {
    double freq_Hz = std::numeric_limits<double>::quiet_NaN();  // Frequency in Hz
    double step_size = 10.0;                                    // Step size for ray tracing [km]
    bool correction = true;                  // Apply correction to the ray tracing
    bool fine_correction = false;            // Apply fine correction to the ray tracing
    double cutoff_r = 4 * RE;                // Cutoff radius for ray tracing [km]
    double gradn_dx = 1.0;                   // Gradient step size for refractive index [km]
    std::string integ_method = "Euler";      // Integration method ("RK4" or "Euler")
    std::string correction_method = "grid";  // Correction method ("grid" or "newton")
    double kp = -1;                          // Kp index for the ionosphere model
    double rz12 = -1.0;                      // IRI R12 (Rz12) sunspot index: >0 uses the value
                                             // (0<R12<=200); -1 historical/projected (+storm
                                             // model); -2 historical/projected (no storm model)
    bool use_fortran_gcpm = true;            // Use Fortran for ray tracing
    double corr_tol = 1.0;                   // Correction tolerance [m] for ray tracing
    bool compute_higher_order = true;        // Compute second-order delays
    bool use_adaptive_step = true;           // Use adaptive step size for ray tracing
    bool straight_ray = false;               // Assume straight ray path
  };

  struct PathProfile {
    // per each timestep
    VecXd s;             // Path length
    VecXd tec_section;   // Total Electron density for each section [TECU]
    VecXd r;             // Geocentric radial distance
    MatXd pos_eci;       // Position in Cartesian coordinates
    VecXd az_dir;        // Azimuth direction
    VecXd el_dir;        // Elevation direction
    VecXd dist_to_line;  // Distance to the straight line

    // final values
    double sf;                   // Final distance of the ray [km]
    Vec3d dir_start;             // Direction vector at the start
    Vec3d dir_end;               // Direction vector at the end
    Vec3d final_pos;             // Final position in Cartesian coordinate
    Vec3d corr_final_pos_err;    // Final position error in Cartesian coordinate due
                                 // to insufficient correction
    double corr_final_time_err;  // Final time error [s] due to insufficient
                                 // correction
    double t_tx;                 // Transmission epoch (seconds since J2000)
    double t_rx;                 // Reception epoch (seconds since J2000)
    double prop_time_total;      // Total propagation time [s] (technically t_rx -
                                 // t_tx, but computed with care in float limits)
    double dist_bend_m;          // Distance bent [m]
    double dist_straight_km;     // Straight distance [km]
    double total_delay_m;        // Total delay [m]
    double tec_delay_m;          // TEC delay [m]
    double tec_delay_bend_m;     // additional TEC due to ray bending [m] (part of
                                 // tec_delay_m)
    double second_delay_m;       // Second order delay [m]
    double third_delay_m;        // Third order delay [m] (if computed)
    double tecu;                 // Total electron content [TECU]
    double tecu_bend;            // TEC bend [TECU]
    double max_sep_line_m;       // Maximum separation from the straight line [m]
  };

  bool is_method_rk4(const std::string& method);

  /**
   * @brief Convert azimuth and elevation angles to a unit vector
   * @param az Azimuth angle in radians
   * @param el Elevation angle in radians
   * @return Unit vector in Cartesian coordinates
   */
  Vec3d azel_to_unitvec(double az, double el);

  /**
   * @brief Convert a unit vector to azimuth and elevation angles
   * @param x Unit vector in Cartesian coordinates
   * @return Azimuth and elevation angles in radians
   */
  Vec2d unitvec_to_azel(const Vec3d& x);

  /**
   * @brief Compute the refractive index of the ionosphere at a given position
   * @param t_j2000 Time in seconds since J2000 epoch
   * @param pos_geo Position in geocentric coordinates (ECEF)
   * @param config Ray trace configuration
   * @param debug Debug flag to print additional information
   * @return Refractive index at the given position
   */
  double compute_ne(double t_j2000, const Vec3d& pos_geo, RayTraceConfig config,
                    bool debug = false);

  /**
   * @brief Compute the ionospheric parameters at a given time in J2000 format
   * @param t_j2000 Time in seconds since J2000 epoch
   * @param kp Kp index for the ionosphere model (default is -1.0
   *  to use the default value)
   * @return Vector containing [f107, rz12, hmf2_km]
   */
  Vec3d get_iono_params(double t_j2000, double kp = -1.0);

  /**
   * @brief Compute the magnetic field strength at a given position
   * @param t_j2000 Time in seconds since J2000 epoch
   * @param pos_geo Position in geocentric coordinates (ECEF)
   * @param config Ray trace configuration
   * @param debug Debug flag to print additional information
   * @return Magnetic field strength vector in Cartesian coordinates [nT]
   */
  Vec3d compute_B(double t_j2000, const Vec3d& pos_geo, RayTraceConfig config, bool debug = false);

  /**
   * @brief Compute the refractive index of the ionosphere at a given position
   * @param ne_m3 Electron density in m^-3
   * @param B Magnetic field strength in nT
   * @param costheta Cosine of the angle between the magnetic field and the ray
   * @param freq_Hz Frequency in Hz
   * @return Refractive index at the given electron density and frequency
   */
  double refractive_index_neB(double ne_m3, double freq_Hz, double B, double costheta,
                              bool compute_higher_order = true);

  /**
   * @brief Compute the refractive index of the ionosphere at a given time and
   * position
   * @param t Time in seconds since J2000 epoch
   * @param x Position in Cartesian coordinates (ECEF)
   * @param shat Direction vector of the ray (normalized)
   * @param config Ray trace configuration
   * @param debug Debug flag to print additional information
   * @return Refractive index at the given time and position
   */
  double refractive_index(double t, const Vec3d& x, const Vec3d& shat, RayTraceConfig config,
                          bool debug = false);

  /**
   * @brief Compute the gradient of the refractive index
   * @param t Time in seconds since J2000 epoch
   * @param x Position in Cartesian coordinates (ECEF)
   * @param config Ray trace configuration
   * @param debug Debug flag to print additional information
   */
  Vec3d gradient_n(double t, const Vec3d& x, RayTraceConfig config, bool debug = false);

  /**
   * @brief Compute the derivative of the ray at a given point
   * @param s Path length at the point
   * @param z State vector containing position, direction, time, and TEC
   * @param t_epoch_tx Transmission epoch in seconds since J2000
   * @param config Ray trace configuration
   * @param store_vec  n (1) and grad_n (3) storage
   * @param debug Debug flag to print additional information
   * @return Derivative of the ray at the given point
   */
  VecXd ray_derivative(double s, const VecXd& z, double t_epoch_tx, RayTraceConfig config,
                       VecXd& store_vec, bool use_precomputed_vals, bool debug = false);

  /**
   * @brief Integrate the ray using a numerical method
   * @param z State vector containing position, direction, time,
   * and TEC
   * @param s Path length at the point
   * @param h Step size for integration
   * @param t_epoch_tx Transmission epoch in seconds since J2000
   * @param config Ray trace configuration
   * @param debug Debug flag to print additional information
   * @return Updated state vector after integration
   */
  void integ_step(VecXd& z, double& s, double h, double t_epoch_tx, RayTraceConfig config,
                  VecXd& store_vec, bool use_precomputed_vals, bool debug = false);

  double adjust_stepsize(double r, const VecXd& store_vec);

  /**
   * @brief Propagate a ray from transmitter to receiver
   * @param t_prop Propagation time in seconds (=t_rx - t_tx)
   * @param dir Initial direction vector of the ray (normalized)
   * @param sf Final distance of the ray (computed by correction)
   * @param x1 Position of the transmitter in Cartesian coordinates
   * @param x2 Position of the receiver in Cartesian coordinates
   * @param t_rx Reception time in seconds since J2000 epoch
   * @param config Ray trace configuration
   * @param azel_diff Azimuth and elevation difference for correction
   * @param store_mat  Matrix to store n, gradn, ne_m3, and B_dot_s for each step
   * @param use_precomputed_densities Use precomputed dzds for faster propagation
   * @return Final state vector containing position, direction,
   * time, and TEC
   */
  VecXd propagate_ray(double t_prop, const Vec3d& dir, double sf, const Vec3d& x1, const Vec3d& x2,
                      double t_rx, RayTraceConfig config, MatXd& store_mat,
                      bool use_precomputed_densities = false);

  /**
   * @brief Trace a ray from transmitter to receiver
   * @param epoch_utc_rx Reception time in seconds since J2000 epoch
   * @param x1 Position of the transmitter in Cartesian coordinates
   * @param x2 Position of the receiver in Cartesian coordinates
   * @param config Ray trace configuration
   * @param debug_prop Debug flag to print additional information for propagation
   * @param debug_corr Debug flag to print additional information for correction
   * @return PathProfile containing the path profile information
   */
  PathProfile trace_ray(double epoch_utc_rx, const Vec3d& x1, const Vec3d& x2,
                        RayTraceConfig config, bool debug_prop = false, bool debug_corr = false);

  /**
   * @brief Compute the path profile for a ray from transmitter to receiver
   * @param t_prop Propagation time in seconds (=t_rx - t_tx)
   * @param dir Inital direction vector of the ray (normalized, computed by
   * correction)
   * @param sf Final distance of the ray (computed by correction)
   * @param txpos Function to get the position of the transmitter at time t_tx
   * @param x2 Position of the receiver in Cartesian coordinates
   * @param t_rx Reception time in seconds since J2000 epoch
   * @param config Ray trace configuration
   * @param debug Debug flag to print additional information
   * @return PathProfile containing the path profile information
   */
  PathProfile compute_path_profile(double t_prop, const Vec3d& dir, double sf, const Vec3d& x1,
                                   const Vec3d& x2, double t_rx, RayTraceConfig config,
                                   MatXd& store_mat, bool use_precomputed_densities = false,
                                   bool debug = false);

  /** * @brief Compute the TEC for a straight line distance
   * @param t_tx Transmission epoch in seconds since J2000 epoch
   * @param initial_pos Initial position of the ray in Cartesian coordinates
   * @param final_pos Final position of the ray in Cartesian coordinates
   * @param config Ray trace configuration
   * @return TEC for the straight line distance in meters
   * */
  double compute_tec_straight(double t_tx, const Vec3d& initial_pos, const Vec3d& final_pos,
                              RayTraceConfig config);

  /** * @brief Compute the TEC for a section of the ray
   * @param i Index of the section
   * @param t_tx Transmission epoch in seconds since J2000 epoch
   * @param initial_pos Initial position of the ray in Cartesian coordinates
   * @param dir Direction vector of the ray (normalized)
   * @param config Ray trace configuration
   * @return TEC for the section in meters
   * */
  double compute_tec_section(int i, double t_tx, const Vec3d& initial_pos, const Vec3d& dir,
                             const RayTraceConfig& config);

  void correct_proptime(double& prop_time, Vec3d& dir, double sf, const Vec3d& x1, const Vec3d& x2,
                        double t_rx, RayTraceConfig config, VecXd& zf);

  void correct_dir_neldermead(double& t_prop, Vec3d& dir, double sf, const Vec3d& x1,
                              const Vec3d& x2, double t_rx, RayTraceConfig config, VecXd& zf,
                              MatXd& store_mat,
                              bool coarse_mode = false,  // Use precomputed dzds
                              bool debug = false);
  /**
   * @brief Correct the ray tracing result using a Newton method
   * @param iter Current iteration number
   * @param t_prop Propagation time in seconds
   * @param dir Initial direction vector of the ray (normalized)
   * @param sf Path length of the ray (computed by correction)
   * @param txpos Function to get the position of the transmitter at time t_tx
   * @param x2 Position of the receiver in Cartesian coordinates
   * @param t_rx Reception time in seconds since J2000 epoch
   * @param config Ray trace configuration
   * @param zf Final state vector containing the ray tracing results
   * @param debug Debug flag to print additional information
   *  */
  void correct_dir_newton(double& t_prop, Vec3d& dir, double sf, const Vec3d& x1, const Vec3d& x2,
                          double t_rx, RayTraceConfig config, VecXd& zf, bool debug = false);

}  // namespace pecsim
