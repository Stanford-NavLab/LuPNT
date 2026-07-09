/**
 * @file raytrace.cpp
 * @author Keidai Iiyama
 * @brief This file contains the implementation of the ray tracing function
 * @version 0.1
 * @date 2025-02-17
 *
 * @copyright Copyright (c) 2025
 *
 */
#include "lupnt/environment/plasma/tec/raytrace.h"

#include <omp.h>

#include "lupnt/environment/plasma/core/progress_bar.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"
#include "lupnt/environment/plasma/gcpm/gcpm_interface.h"
#include "lupnt/environment/plasma/igrf/igrf_interface.h"
#include "lupnt/environment/plasma/tec/neldermead.h"

namespace pecsim {

  static const int STATE_SIZE = 10;        // Size of the state vector
  static const int STORE_VEC_SIZE = 11;    // Size of the storage vector
  static const double p_coeff = 40.3;      // Coefficient for first-order delay
  static const double q_coeff = 2.2566e3;  // For second-order delay (for B in nT)
  static const double u1_coeff = 2437;
  static const double u2_coeff = 4.74e4;  // For third-order delay (for B in nT)

  bool is_method_rk4(const std::string& method) { return (method == "RK4" || method == "rk4"); }

  Vec3d azel_to_unitvec(double az_rad, double el_rad) {
    // Convert azimuth and elevation angles to a unit vector
    double cos_el = std::cos(el_rad);
    Vec3d uv;
    uv[0] = cos_el * std::cos(az_rad);
    uv[1] = cos_el * std::sin(az_rad);
    uv[2] = std::sin(el_rad);
    return uv.normalized();  // Ensure the vector is normalized
  }

  Vec2d unitvec_to_azel(const Vec3d& x) {
    // Convert a unit vector to azimuth and elevation angles
    Vec3d xn = x / x.norm();  // Ensure x is a unit vector
    double az_rad = std::atan2(xn(1), xn(0));
    double el_rad = std::asin(xn(2));
    return Vec2d(az_rad, el_rad);
  }

  double compute_ne(double t_j2000, const Vec3d& pos_geo, RayTraceConfig config, bool debug) {
    // Extract configuration parameters
    double kp = config.kp;  // Kp index for the ionosphere model
    if (kp < 0) {
      throw std::runtime_error("Kp index is not set. Please set the Kp index.");
    }

    // First check if iri model is set
    if (get_iri_model() == IRIModel::NONE) {
      throw std::runtime_error("IRI model is not set. Please set the IRI model.");
    }

    double alatr = 0.0;
    double along = 0.0;
    double r = 0.0;
    std::array<double, 3> pos_gei_arr = {pos_geo[0] / RE, pos_geo[1] / RE, pos_geo[2] / RE};
    std::array<double, 3> pos_sm;

    // convert to magnetic coordinates
    double mjd = tj2000_to_mjd(t_j2000);
    DateTime datetime = mjd_to_datetime(mjd);
    std::array<int, 2> itime = datetime_to_itime(datetime);

    geo_to_sm(itime, pos_gei_arr, pos_sm);
    cart_to_pol(pos_sm, alatr, along, r);
    double amlt = 12.0 + along / AMLTRAD;
    if (amlt > 24.0) amlt -= 24.0;

    std::vector<double> out;

    if (config.use_fortran_gcpm)
      out = gcpm_v24_fortran(datetime, r, amlt, alatr, kp);
    else
      out = gcpm_v24(datetime, r, amlt, alatr, kp);

    double ne = out[0];  // the first element is the electron density

    if (debug) {
      std::cout << "    [compute_ne] " << std::endl;
      std::cout << "      Time (MJD): " << mjd << " itime: " << itime[0] << ", " << itime[1]
                << "  Datetime: " << datetime.year << ", " << datetime.doy << ", " << datetime.hour
                << ", " << datetime.min << ", " << datetime.sec << std::endl;
      std::cout << "      ALAT: " << alatr * RAD2DEG << " deg"
                << "  ALONG: " << along * RAD2DEG << " deg"
                << "  AMLT: " << amlt << std::endl;
      std::cout << "      ne (cm^-3): " << ne << std::endl;
    }

    return ne;
  }

  Vec3d get_iono_params(double t_j2000, double kp) {
    // Get the ionospheric parameters at a given time
    double mjd = tj2000_to_mjd(t_j2000);
    DateTime datetime = mjd_to_datetime(mjd);

    // dummy variables for the function call
    double r = 1.05;
    double alatr = 0.0;  // Geomagnetic latitude in radians
    double amlt = 23.0;  // Geomagnetic local time in hours

    // Call the GCPM function to get the ionospheric parameters
    std::vector<double> out = gcpm_v24(datetime, r, amlt, alatr, kp);

    Vec3d iono_params;
    iono_params[0] = out[4];  // f107 index
    iono_params[1] = out[5];  // rz12 index
    iono_params[2] = out[6];  // hmf2_km (F2 peak height in km)

    return iono_params;
  }

  Vec3d compute_B(double t_j2000, const Vec3d& pos_geo, RayTraceConfig config, bool debug) {
    // Compute the magnetic field strength at a given position
    // This is a placeholder function; actual implementation may vary
    // depending on the magnetic field model used.
    // For now, we return a constant value for demonstration purposes.
    double lat_rad, lon_rad, r;
    std::array<double, 3> pos_geo_arr = {pos_geo[0], pos_geo[1], pos_geo[2]};

    double mjd = tj2000_to_mjd(t_j2000);
    DateTime datetime = mjd_to_datetime(mjd);
    double decimal_year = datetime.year + (datetime.doy - 1) / 365.25;

    cart_to_pol(pos_geo_arr, lat_rad, lon_rad, r);
    std::vector<double> B = igrf14(lat_rad, lon_rad, r, decimal_year);

    Vec3d B_vec = {
        B[0],  // Bx
        B[1],  // By
        B[2]   // Bz
    };

    if (debug) {
      std::cout << "    [compute_B] " << std::endl;
      std::cout << "      Time (MJD): " << mjd << "  Datetime: " << datetime.year << ", "
                << datetime.doy << ", " << datetime.hour << ", " << datetime.min << ", "
                << datetime.sec << std::endl;
      std::cout << "      Position (Geo): (" << pos_geo[0] << ", " << pos_geo[1] << ", "
                << pos_geo[2] << ") km" << std::endl;
      std::cout << "      Magnetic Field (nT): (" << B_vec[0] * 1e9 << ", " << B_vec[1] * 1e9
                << ", " << B_vec[2] * 1e9 << ") nT" << std::endl;
    }

    return B_vec;
  }

  double refractive_index_neB(double ne_m3, double freq_Hz, double B, double cos_theta,
                              bool compute_q) {
    // Calculate the refractive index based on electron density and frequency
    // ne_m3: electron density in m^-3
    // freq_Hz: frequency in Hz
    double f2 = freq_Hz * freq_Hz;
    double f3 = f2 * freq_Hz;
    double f4 = f3 * freq_Hz;
    double p = p_coeff * ne_m3;
    double q = 0.0;  // second order correction term
    double u = 0.0;  // third order
    if (compute_q) {
      // Calculate the correction term q based on the frequency
      // This is a simplified model; actual implementation may vary
      q = q_coeff * ne_m3 * B * cos_theta;  // B is in nT, ne_m3 in m^-3
      u = u1_coeff * ne_m3 * ne_m3 + u2_coeff * ne_m3 * B * B * (1 + cos_theta * cos_theta);
    }

    double n = 1 - (p / f2 + q / (2 * f3) + u / (3 * f4));  // phase refractive index

    return n;  // Return the refractive index
  }

  double refractive_index(double t, const Vec3d& x, const Vec3d& shat, RayTraceConfig config,
                          bool debug) {
    // Extract configuration parameters
    double freq_Hz = config.freq_Hz;  // Frequency in Hz

    double ne_cm3 = compute_ne(t, x, config, debug);
    double ne_m3 = ne_cm3 * 1e6;  // Convert from cm^-3 to m^-3
    double B = 0.0, cos_theta = 0.0;

    bool compute_higher_order = config.compute_higher_order;
    if (compute_higher_order) {
      // Compute the magnetic field at the position x
      Vec3d Bvec = compute_B(t, x, config, debug);
      double Bdots = Bvec.dot(shat);
      B = Bvec.norm();
      cos_theta = Bdots / B;
    }

    double n = refractive_index_neB(ne_m3, freq_Hz, B, cos_theta, config.compute_higher_order);

    return n;
  }

  Vec3d gradient_n(double t, const Vec3d& x, RayTraceConfig config, bool debug) {
    Vec3d gradn;
    double dx = config.gradn_dx;  // Position deviation [km] to compute the
                                  // derivative of the ion distribution

    // Extract configuration parameters
    double freq_Hz = config.freq_Hz;  // Frequency in Hz

    if (debug) {
      std::cout << "    [gradient_n]  t: " << t << std::endl;
    }

    for (int i = 0; i < 3; ++i) {
      Vec3d x_forward = x;
      Vec3d x_backward = x;
      x_forward(i) += dx;
      x_backward(i) -= dx;

      Vec3d shat = Vec3d::Zero();  // just use dummy

      double nf = refractive_index(t, x_forward, shat, config, false);
      double nb = refractive_index(t, x_backward, shat, config, false);
      gradn(i) = (nf - nb) / (2 * dx);

      if (debug) {
        std::vector<std::string> xyz_str = {"x", "y", "z"};
        std::cout << "       " << xyz_str[i] << "|  grad: " << gradn(i) << std::endl;
      }
    }
    return gradn;
  }

  VecXd ray_derivative(double s, const VecXd& z, double t_epoch_tx, RayTraceConfig config,
                       VecXd& store_vec, bool use_precomputed_vals, bool debug) {
    // Extract configuration parameters
    double freq_Hz = config.freq_Hz;    // Frequency in Hz
    double dx = config.gradn_dx;        // Position deviation [km] to compute the
                                        // derivative of the ion distribution
    double cutoff_r = config.cutoff_r;  // Cutoff radius for ray tracing [km]

    // States
    Vec3d x = z.segment<3>(0);
    Vec3d s_hat = z.segment<3>(3);
    if (s_hat.norm() > 0) s_hat.normalize();
    double t_epoch = z(6) + t_epoch_tx;
    double r = x.norm();

    if (store_vec.size() != STORE_VEC_SIZE) {
      if (use_precomputed_vals) {
        throw std::runtime_error("Store vector size mismatch. Expected: "
                                 + std::to_string(STORE_VEC_SIZE)
                                 + ", but got: " + std::to_string(store_vec.size()));
      } else {
        store_vec.resize(STORE_VEC_SIZE);
      }
    }

    // flags
    bool compute_ho = config.compute_higher_order;

    // refractive index n and gradn computetation
    double n = 1.0;
    Vec3d grad_n;
    double ne_m3 = 0.0, B = 0.0, cos_theta = 0.0;

    if (r >= cutoff_r) {
      // skip computation if radius is over cutoff threshold
      grad_n.setZero();
      ne_m3 = 0.0;
      B = 0.0;
      cos_theta = 0.0;

    } else {
      if (use_precomputed_vals) {
        n = store_vec(0);
        grad_n = store_vec.segment<3>(1);
        ne_m3 = store_vec(4);
        B = store_vec(5);
        cos_theta = store_vec(6);

      } else {
        double ne_cm3 = compute_ne(t_epoch, x, config, debug);
        ne_m3 = ne_cm3 * 1e6;  // Convert from cm^-3 to m^-3

        if (compute_ho) {
          Vec3d B_vec = compute_B(t_epoch, x, config, debug);  // in nT
          double B_dot_s = B_vec.dot(s_hat);                   // in nT m^3
          B = B_vec.norm();                                    // in nT
          cos_theta = B_dot_s / B;                             // cos(theta
        } else {
          B = 0.0;
          cos_theta = 0.0;
        }

        n = refractive_index_neB(ne_m3, freq_Hz, B, cos_theta, compute_ho);
        if (config.straight_ray) {
          grad_n.setZero();
        } else {
          grad_n = gradient_n(t_epoch, x, config, debug);
        }
      }
    }

    Vec3d dshat_ds = Vec3d::Zero();

    if (config.straight_ray) {
      dshat_ds.setZero();
    } else {
      dshat_ds = (grad_n - s_hat * grad_n.dot(s_hat)) / n;
    }

    // Compute the derivatives -----------------------------------
    VecXd dz(STATE_SIZE);
    dz.segment<3>(0) = s_hat;
    dz.segment<3>(3) = dshat_ds;

    // We only track the time increase due to distance here due
    // to machine precision in tracking TEC delays (=p_coeff *
    // ne_m3 / (freq_Hz  * freq_Hz) / C) = ne_m3 * p_coeff * 1e-23 )
    // -> We do not want to keep epochs and TEC delays in the same vector
    dz(6) = 1.0 / C;
    dz(7) = ne_m3 * 1000;  // since the s is in km, to get TEC in m^-2
    dz(8) = q_coeff * ne_m3 * B * cos_theta
            * 1000;  // For the second order delay (since the s is in km)
    dz(9) = 1000
            * (u1_coeff * ne_m3 * ne_m3 + u2_coeff * ne_m3 * B * B * (1 + cos_theta * cos_theta));

    // storage ------------------------------------------------------
    store_vec(0) = n;
    store_vec.segment<3>(1) = grad_n;
    store_vec(4) = ne_m3;
    store_vec(5) = B;
    store_vec(6) = cos_theta;
    // Store position and time for paralllel computation of n, gradn
    store_vec.segment<3>(7) = x;
    store_vec(10) = s;

    if (debug) {
      std::cout << "    [ray_derivative] dz = " << dz.transpose() << std::endl;
    }

    return dz;
  }

  void integ_step(VecXd& z, double& s, double h, double t_epoch_tx, RayTraceConfig config,
                  VecXd& store_vec, bool use_precomputed_vals, bool debug) {
    // unpack
    std::string method = config.integ_method;

    if (debug) {
      std::cout << "  [integ_step]" << std::endl;
    }

    int store_size = STORE_VEC_SIZE;
    if (is_method_rk4(method)) {
      store_size *= 4;  // For RK4, we need to store 4
    }

    if (store_vec.size() != store_size) {
      if (use_precomputed_vals) {
        throw std::runtime_error("Store vector size mismatch. Expected: "
                                 + std::to_string(store_size)
                                 + ", but got: " + std::to_string(store_vec.size()));
      } else {
        // Initialize the store vector if not using precomputed values
        store_vec.resize(store_size);
        store_vec.setZero();  // Initialize to zero
      }
    }

    if (is_method_rk4(method)) {
      VecXd store_vec1(STORE_VEC_SIZE), store_vec2(STORE_VEC_SIZE), store_vec3(STORE_VEC_SIZE),
          store_vec4(STORE_VEC_SIZE);
      if (use_precomputed_vals) {
        store_vec1 = store_vec.segment<STORE_VEC_SIZE>(0);
        store_vec2 = store_vec.segment<STORE_VEC_SIZE>(STORE_VEC_SIZE);
        store_vec3 = store_vec.segment<STORE_VEC_SIZE>(2 * STORE_VEC_SIZE);
        store_vec4 = store_vec.segment<STORE_VEC_SIZE>(3 * STORE_VEC_SIZE);
        h = store_vec4(10) - store_vec1(10);  // Use the last step size
      } else {
        store_vec1.setZero();
        store_vec2.setZero();
        store_vec3.setZero();
        store_vec4.setZero();
      }

      VecXd k1 = ray_derivative(s, z, t_epoch_tx, config, store_vec1, use_precomputed_vals, debug);
      VecXd z2 = z + 0.5 * h * k1;

      VecXd k2 = ray_derivative(s + 0.5 * h, z2, t_epoch_tx, config, store_vec2,
                                use_precomputed_vals, false);
      z2 = z + 0.5 * h * k2;

      VecXd k3 = ray_derivative(s + 0.5 * h, z2, t_epoch_tx, config, store_vec3,
                                use_precomputed_vals, false);
      z2 = z + h * k3;

      VecXd k4
          = ray_derivative(s + h, z2, t_epoch_tx, config, store_vec4, use_precomputed_vals, false);
      z += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

      // Store the derivatives
      store_vec.segment<STORE_VEC_SIZE>(0) = store_vec1;
      store_vec.segment<STORE_VEC_SIZE>(STORE_VEC_SIZE) = store_vec2;
      store_vec.segment<STORE_VEC_SIZE>(2 * STORE_VEC_SIZE) = store_vec3;
      store_vec.segment<STORE_VEC_SIZE>(3 * STORE_VEC_SIZE) = store_vec4;

    } else if ((method == "EULER") || (method == "euler") || (method == "Euler")) {
      // Save the stored step size BEFORE ray_derivative overwrites slot 10
      // with the current arc length. Slot 10 holds step_size in non-precomputed
      // mode (written below) and is used to replay the same step sizes during
      // precomputed coarse-NM correction.
      if (use_precomputed_vals) {
        h = store_vec(10);  // save stored step_size before ray_derivative clobbers it
      }
      VecXd dz = ray_derivative(s, z, t_epoch_tx, config, store_vec, use_precomputed_vals, debug);
      if (!use_precomputed_vals) {
        store_vec(10) = h;  // Store the step size for future precomputed replay
      }
      z += h * dz;
    } else {
      throw std::runtime_error("Unknown integration method: " + method);
    }

    // Update the path length
    s += h;  // Increment the path length
  }

  double adjust_stepsize(double r, const VecXd& store_vec) {
    // Adjust the step size based on the radial distance and the stored values
    // This is a placeholder function; actual implementation may vary
    // depending on the specific requirements of the ray tracing algorithm.
    double step_size = 20.0;  // Default step size in km
    double alt = r - RE;      // Altitude in km

    if (alt < 1000) {
      step_size = 10.0;  // Smaller step size for closer distances
    } else if (alt < 4000) {
      step_size = 20.0;  // Default step size for medium distances
    } else {
      step_size = 100.0;  // Maximum step size for very far distances
    }

    return step_size;
  }

  VecXd propagate_ray(double t_prop, const Vec3d& dir, double sf, const Vec3d& x1, const Vec3d& x2,
                      double t_rx, RayTraceConfig config, MatXd& store_mat,
                      bool use_precomputed_vals) {
    // Extract configuration parameters
    double freq_Hz = config.freq_Hz;                 // Frequency in Hz
    double step_size = config.step_size;             // Step size for ray tracing [km]
    bool correction = config.correction;             // Apply correction to the ray tracing
    double cutoff_r = config.cutoff_r;               // Cutoff radius for ray tracing [km]
    std::string integ_method = config.integ_method;  // Integration method ("RK4" or "Euler")
    double kp = config.kp;                           // Kp index for the ionosphere model

    double s = 0.0;
    double t_tx = t_rx - t_prop;  // Transmission time in seconds since J2000

    // initialize state
    VecXd z(STATE_SIZE);
    z.segment<3>(0) = x1;
    z.segment<3>(3) = dir;
    z(6) = 0.0;  // Time since transittion epoch
    z(7) = 0.0;  // TEC
    z(8) = 0.0;  // Second order delay (if computed)
    z(9) = 0.0;  // Third order delay (if computed)

    // When adaptive stepping is active, step_size can drop to 10 km (the
    // minimum returned by adjust_stepsize) even if config.step_size is large.
    // Allocate for worst-case step count so sidx never overruns store_mat.
    const double min_adaptive_step_p = 10.0;
    double alloc_step_p = (config.use_adaptive_step && !use_precomputed_vals)
                              ? std::min(step_size, min_adaptive_step_p)
                              : step_size;
    int n_steps = static_cast<int>(sf / alloc_step_p) + 2;
    int sidx = 0;               // Step index for debugging output
    bool cutoff_hit_p = false;  // Set true when the ray exits via the cutoff radius
    Vec3d pos_now, dir_now;
    double r = 0.0;           // Geocentric radial distance
    double r_prev = 50 * RE;  // Previous radial distance for comparison (set to
                              // a large value to ensure the first step is
                              // always executed)
    // update storage
    // n(1), grad_n(3), ne_m3(1), B_dot_s(1)
    if (!use_precomputed_vals) {
      if (is_method_rk4(integ_method)) {
        // If using RK4, we need to store the state vector for each step
        store_mat.resize(n_steps, STORE_VEC_SIZE * 4);
      } else {
        store_mat.resize(n_steps, STORE_VEC_SIZE);
      }
    }

    while (s < sf) {
      if ((sf - s) < step_size) {
        step_size = sf - s;
      }

      VecXd z_orig = z;  // Store the original state vector for debugging

      // Compute the derivative and store it in store_mat
      VecXd store_vec;
      if (use_precomputed_vals) {
        if (sidx >= store_mat.rows()) break;  // perturbed ray is longer than precomputed path
        store_vec = store_mat.row(sidx).transpose();
      } else {
        if (is_method_rk4(integ_method)) {
          // If using RK4, we need to store the state vector for each step
          store_vec.resize(STORE_VEC_SIZE * 4);
        } else {
          store_vec.resize(STORE_VEC_SIZE);
        }
      }

      // Integrate the ray using the specified method
      integ_step(z, s, step_size, t_tx, config, store_vec, use_precomputed_vals, false);

      // Store
      store_mat.row(sidx) = store_vec.transpose();

      // Normalize the direction vector
      z.segment<3>(3) = z.segment<3>(3).normalized();

      // Compute variables
      pos_now = z.segment<3>(0);
      r = pos_now.norm();
      dir_now = z.segment<3>(3);

      // Decide next step based on the radial distance
      if ((config.use_adaptive_step) && (!use_precomputed_vals)) {
        step_size = adjust_stepsize(r, store_vec);
      }

      // If it exits the cutoff radius, stop the propagation
      if ((r >= cutoff_r) && (r > r_prev)) {
        // Shoot for the final position ------------------------------
        // Find sf that minimize |fpos - x2|, where fpos = pos_now + dir_now *
        // (sf - s)
        Vec3d delta_x2 = x2 - pos_now;
        double delta_s = delta_x2.dot(dir_now);
        sf = s + delta_s;
        Vec3d final_pos = pos_now + dir_now * (sf - s);
        double computed_prop_time = z(6) + (sf - s) / C;  // Final time in seconds since J2000

        z.segment<3>(0) = final_pos;  // Update the position in the state vector
        z(6) = computed_prop_time;    // Update the time in the state vector

        cutoff_hit_p = true;
        break;  // Exit the loop after reaching the cutoff radius
      } else {
        sidx++;
      }

      r_prev = r;  // Store the previous radial distance
    }

    // Trim store_mat to the rows that were actually written so that subsequent
    // precomputed calls (Nelder-Mead) cannot read uninitialized rows (h=0 would
    // cause an infinite loop since s never advances, and sidx would overflow).
    // Cutoff exit: last write was at row sidx (sidx++ was NOT called before break).
    // Non-cutoff exit: sidx++ ran after each write, so last write was at row sidx-1.
    if (!use_precomputed_vals) {
      int valid_rows = cutoff_hit_p ? (sidx + 1) : sidx;
      int store_cols = is_method_rk4(integ_method) ? STORE_VEC_SIZE * 4 : STORE_VEC_SIZE;
      store_mat.conservativeResize(valid_rows, store_cols);
    }

    return z;
  }

  double compute_tec_straight(double t_tx, const Vec3d& initial_pos, const Vec3d& final_pos,
                              RayTraceConfig config) {
    double dist = (final_pos - initial_pos).norm();  // Distance in km
    double t_prop = dist / C;                        // Propagation time in seconds
    int n_steps = static_cast<int>(dist / config.step_size) + 1;
    int sidx = 0;                                        // Step index for debugging output
    VecXd dir = (final_pos - initial_pos).normalized();  // Direction vector
    VecXd tec_section = VecXd::Zero(n_steps);            // TEC sections

    // First, compute the time and position vectors alogn the ray path
    for (int i = 0; i < n_steps; ++i) {
      tec_section(i) = compute_tec_section(i, t_tx, initial_pos, dir, config);
    }

    // Sum the TEC sections to get the total TEC
    double tec_total = tec_section.sum();

    return tec_total;  // Return the total TEC in m^-2
  }

  double compute_tec_section(int i, double t_tx, const Vec3d& initial_pos, const Vec3d& dir,
                             const RayTraceConfig& config) {
    double step_size = config.step_size;  // Step size for ray tracing [km]
    double s = i * step_size;             // Path length in km
    double t = t_tx + s / C;              // Time in seconds since J2000

    // Compute the position along the ray path
    Vec3d pos_now = initial_pos + dir * s;

    if (pos_now.norm() > config.cutoff_r) {
      // If the position is outside the cutoff radius, return 0
      return 0.0;
    }

    // Compute the electron density at this position
    double ne_cm3 = compute_ne(t, pos_now, config, false);
    double ne_m3 = ne_cm3 * 1e6;  // Convert from cm^-3 to m^-3

    // Compute the TEC contribution for this step
    double tec_section = ne_m3 * 1000 * step_size;  // in m^-2

    return tec_section;
  }

  PathProfile compute_path_profile(double t_prop, const Vec3d& dir, double sf, const Vec3d& x1,
                                   const Vec3d& x2, double t_rx, RayTraceConfig config,
                                   MatXd& store_mat, bool use_precomputed_densities, bool debug) {
    // Extract configuration parameters
    double freq_Hz = config.freq_Hz;
    double step_size = config.step_size;             // Step size for ray tracing [km]
    bool correction = config.correction;             // Apply correction to the ray tracing
    double cutoff_r = config.cutoff_r;               // Cutoff radius for ray tracing [km]
    std::string integ_method = config.integ_method;  // Integration method ("RK4" or "Euler")
    double kp = config.kp;                           // Kp index for the ionosphere model

    // propagate the ray using RK4
    double s = 0.0;

    // Initial position
    double t_tx = t_rx - t_prop;  // Transmission time in seconds since J2000
    Vec3d initial_pos = x1;

    // initialize state
    VecXd z(STATE_SIZE);
    z.segment<3>(0) = initial_pos;
    z.segment<3>(3) = dir.normalized();
    z(6) = 0.0;  // Time since transmission epoch
    z(7) = 0.0;  // TEC
    z(8) = 0.0;  // Second order delay (if computed)
    z(9) = 0.0;  // Third order delay (if computed)

    // Allocate enough rows for all steps that will actually be taken.
    // Non-precomputed: adaptive stepping can use 10 km steps even if step_size=100,
    //   so allocate based on the minimum adaptive step.
    // Precomputed: step sizes come from the stored h values in store_mat (same
    //   adaptive steps as the original run), NOT from config.step_size.
    //   Using sf/step_size would severely underestimate the row count and cause
    //   OOB writes to pp.s, pp.r, etc., corrupting the heap.
    int n_steps;
    if (use_precomputed_densities) {
      n_steps = static_cast<int>(store_mat.rows()) + 2;
    } else {
      const double min_adaptive_step = 10.0;
      double alloc_step
          = config.use_adaptive_step ? std::min(step_size, min_adaptive_step) : step_size;
      n_steps = static_cast<int>(sf / alloc_step) + 2;
    }

    PathProfile pp;
    pp.s = VecXd::Zero(n_steps);
    pp.r = VecXd::Zero(n_steps);
    pp.tec_section = VecXd::Zero(n_steps);
    pp.az_dir = VecXd::Zero(n_steps);
    pp.el_dir = VecXd::Zero(n_steps);
    pp.pos_eci = MatXd::Zero(n_steps, 3);
    pp.dist_to_line = VecXd::Zero(n_steps);
    pp.dir_start = dir;

    bool cutoff_hit = false;        // Set to true when the ray exits the cutoff radius
    bool debug_integ = false;       // Debug flag for RK4 steps
    if (debug) debug_integ = true;  // Enable detailed debug output for the first step
    int sidx = 0;

    double r = 0.0;          // Geocentric radial distance
    Vec3d pos_now;           // Current position in Cartesian coordinates
    Vec3d dir_now;           // Direction vector
    double tecu = 0.0;       // Total electron content in m^-2
    double tecu_prev = 0.0;  // Previous TEC value for debugging

    Vec3d final_pos = x2;     // Final position; initialized to x2, updated during propagation
    double final_time;        // Final time in seconds since J2000
    double r_prev = 50 * RE;  // Previous radial distance for comparison (set to
                              // a large value to ensure the first step is
                              // always executed)

    // dummy storage vector for ray_derivative
    if (!use_precomputed_densities) {
      if (is_method_rk4(integ_method)) {
        // If using RK4, we need to store the state vector for each step
        store_mat.resize(n_steps, STORE_VEC_SIZE * 4);
      } else {
        store_mat.resize(n_steps, STORE_VEC_SIZE);
      }
    }

    // std::cout << "[compute_path_profile] Starting ray tracing..." << std::endl;
    // ProgressBar pg;
    // pg.start(n_steps);

    while (s < sf) {
      if ((sf - s) < step_size) {
        step_size = sf - s;
      }

      // pg.update(sidx);

      if (debug) {
        if (sidx % 50 == 0) {
          std::cout << " " << std::endl;
          std::cout << "[compute_path_profile] step: " << s << " / " << sf << " (" << sidx << "/"
                    << n_steps << ")" << std::endl;
          std::cout << "  Pos  : " << z.segment<3>(0).transpose() << std::endl;
          std::cout << "  Dir  : " << z.segment<3>(3).transpose() << std::endl;
          std::cout << "  R    : " << z.segment<3>(0).norm() / RE << " RE" << std::endl;
          std::cout << "  Time : " << z(6) << std::endl;
          std::cout << "  TECU: " << z(7) / TECU << std::endl;
          debug_integ = true;  // Enable detailed debug output
        } else {
          debug_integ = false;  // Disable detailed debug output
        }
      }

      // Integration
      if (use_precomputed_densities && sidx >= store_mat.rows()) break;
      VecXd store_vec = store_mat.row(sidx).transpose();
      integ_step(z, s, step_size, t_tx, config, store_vec, use_precomputed_densities, debug_integ);

      // Compute variables
      pos_now = z.segment<3>(0);
      r = pos_now.norm();
      dir_now = z.segment<3>(3).normalized();
      tecu = z(7) / TECU;  // Convert TEC from m^-2 to TECU

      // Decide next step based on the radial distance
      if ((config.use_adaptive_step) && (!use_precomputed_densities)) {
        step_size = adjust_stepsize(r, store_vec);
      }

      Vec2d azel_now = unitvec_to_azel(dir_now);

      // store variables
      pp.s(sidx) = s;
      pp.r(sidx) = r;
      pp.pos_eci.row(sidx) = pos_now.transpose();
      pp.tec_section(sidx) = tecu - tecu_prev;  // TEC section in TECU
      pp.az_dir(sidx) = azel_now(0);            // Azimuth in radians
      pp.el_dir(sidx) = azel_now(1);            // Elevation in radians
      tecu_prev = tecu;

      // If it exits the cutoff radius, stop the propagation
      if ((r >= cutoff_r) && (r > r_prev)) {
        // Shoot for the final position ------------------------------
        // Find sf that minnimize |fpos - x2|, where fpos = pos_now + dir_now *
        // (sf - s)
        Vec3d delta_x2 = x2 - pos_now;
        double delta_s = delta_x2.dot(dir_now);
        double sf_old = sf;  // Store the old sf for debugging
        double sf_new = s + delta_s;
        double ds = sf_new - sf_old;

        sf = s + delta_s;
        final_pos = pos_now + dir_now * (sf - s);
        double prop_time_computed = z(6) + (sf - s) / C;  // propagation time without TEC delay

        z.segment<3>(0) = final_pos;  // Update the position in the state vector
        z(6) = prop_time_computed;    // Update the time in the state vector

        if (debug) {
          std::cout << " " << std::endl;
          std::cout << "[compute_path_profile] Exiting cutoff radius: " << r << " >= " << cutoff_r
                    << std::endl;
          std::cout << " s correction: " << ds << " km" << std::endl;
        }

        // resize the storage vectors to the current size
        cutoff_hit = true;
        pp.s.conservativeResize(sidx + 1);
        pp.r.conservativeResize(sidx + 1);
        pp.tec_section.conservativeResize(sidx + 1);
        pp.pos_eci.conservativeResize(sidx + 1, 3);
        pp.az_dir.conservativeResize(sidx + 1);
        pp.el_dir.conservativeResize(sidx + 1);
        pp.dist_to_line.conservativeResize(sidx + 1);

        break;

      } else {
        // Update states
        sidx++;
      }

      r_prev = r;  // Store the previous radial distance
    }

    // Non-cutoff exit: the loop ended because s >= sf.  sidx was incremented
    // after each write, so entries 0..sidx-1 are valid.  The pp members were
    // over-allocated (worst-case adaptive step count), so trim them now.
    if (!cutoff_hit) {
      pp.s.conservativeResize(sidx);
      pp.r.conservativeResize(sidx);
      pp.tec_section.conservativeResize(sidx);
      pp.pos_eci.conservativeResize(sidx, 3);
      pp.az_dir.conservativeResize(sidx);
      pp.el_dir.conservativeResize(sidx);
      pp.dist_to_line.conservativeResize(sidx);
    }

    // compute tec for straight line distance
    // std::cout << "[compute_path_profile] Computing TEC for straight line..."
    //           << std::endl;

    double tec_straight = z(7);
    if (config.straight_ray == false) {
      tec_straight = compute_tec_straight(t_tx, initial_pos, final_pos, config);
    }

    // convert the distance to straight line
    Vec3d pos_i, w;
    Vec3d v = final_pos - initial_pos;
    VecXd dist_to_line(pp.s.size());
    double max_sep_line = 0.0;  // Maximum distance to the line for debugging
    for (int i = 0; i < pp.s.size(); ++i) {
      pos_i = pp.pos_eci.row(i).transpose();
      w = pos_i - initial_pos;
      dist_to_line(i) = (v.cross(w)).norm() / v.norm();  // Distance to the straight line
      if (dist_to_line(i) > max_sep_line) {
        max_sep_line = dist_to_line(i);  // Update the maximum distance
      }
    }
    pp.dist_to_line = dist_to_line;  // Store the distance to the line

    double freq2 = freq_Hz * freq_Hz;
    double freq3 = freq2 * freq_Hz;
    double freq4 = freq3 * freq_Hz;

    final_pos = z.segment<3>(0);
    double dist = (final_pos - initial_pos).norm();
    double tec_delay_m = p_coeff * z(7) / freq2;  // TEC delay in meters
    double tec_delay_straight_m
        = p_coeff * tec_straight / freq2;    // Straight line TEC delay in meters
    double second_delay_m = z(8) / (freq3);  // Second order delay in meters
    double third_delay_m = z(9) / (freq4);   // Third order delay in meters

    double C_ms = C * 1000;  // Speed of light in m/s

    // Total propagaion time
    double prop_time_total = z(6) + (tec_delay_m + second_delay_m + third_delay_m) / C_ms;

    // Store results
    pp.dir_end = z.segment<3>(3).normalized();
    pp.sf = sf;
    pp.corr_final_pos_err = final_pos - x2;
    pp.corr_final_time_err = prop_time_total - t_prop;
    pp.t_tx = t_tx;
    pp.t_rx = t_rx;                        // Reception time in seconds since J2000
    pp.prop_time_total = prop_time_total;  // Total propagation time in seconds

    pp.dist_straight_km = dist;
    pp.dist_bend_m = (sf - dist) * 1000;      // sf is the total distance traveled
    pp.max_sep_line_m = max_sep_line * 1000;  // Maximum distance to the straight line in meters

    pp.tec_delay_m = tec_delay_m;                              // TEC delay in meters
    pp.tec_delay_bend_m = tec_delay_m - tec_delay_straight_m;  // Distance bend in meters

    pp.second_delay_m = second_delay_m;  // Second order delay in meters
    pp.third_delay_m = third_delay_m;    // Third order delay in meters (if computed)
    pp.total_delay_m = pp.dist_bend_m + pp.tec_delay_m + pp.second_delay_m
                       + pp.third_delay_m;  // Total delay in m

    pp.tecu = z(7) / TECU;  // the last element is the TEC
    pp.tecu_bend = (z(7) - tec_straight) / TECU;

    return pp;
  }

  PathProfile trace_ray(double epoch_utc_rx, const Vec3d& x1, const Vec3d& x2,
                        RayTraceConfig config, bool debug_prop, bool debug_corr) {
    // Unpack the configuration
    double freq_Hz = config.freq_Hz;                 // Frequency in Hz
    double step_size = config.step_size;             // Step size for the ray tracing in km
    bool correction = config.correction;             // Apply correction to the ray tracing
    double cutoff_r = config.cutoff_r;               // Cutoff radius for the ray tracing
    std::string integ_method = config.integ_method;  // Integration method for ray tracing
    double kp = config.kp;  // Kp index for the ionosphere model (default is 0.0)

    // Apply the configured R12 (Rz12) sunspot index to the active IRI model. rz12 < 0 keeps
    // IRI's historical/projected R12; a positive value (0 < R12 <= 200) is used directly.
    if (get_iri_model() == IRIModel::IRI_2007) {
      IRI2007Option opt;
      opt.R12 = config.rz12;
      opt.update_jf_2007();
      set_iri2007_option(opt);
    } else if (get_iri_model() == IRIModel::IRI_2020) {
      IRI2020Option opt;
      opt.R12 = config.rz12;
      opt.update_jf_2020();
      set_iri2020_option(opt);
    }

    // Compute direction and distance
    Vec3d dir = (x2 - x1).normalized();
    double L = (x2 - x1).norm();
    double prop_time = L / C;  // Light time in seconds

    double sf = L;  // Final point of the ray
    VecXd zf(9);    // Final state vector

    std::string correction_method = config.correction_method;  // Correction method
    if (correction_method != "neldermead" && correction_method != "newton") {
      throw std::runtime_error("Unknown correction method: " + correction_method
                               + " - Supported methods are 'naive' and 'newton'.");
    }

    MatXd store_mat;  // Matrix to store the derivative of the state
                      // derivative vector

    if (correction) {
      if (debug_corr) {
        std::cout << "[trace ray] Start Coarse Update -------------------------" << std::endl;
      }

      // First, run several rounds of coarse correction, fixing the derivative
      // of the direction vector to avoid computing refractive index at each
      // step
      bool coarse_mode = true;                  // Coarse mode for the first correction
      double pos_err = 1e6;                     // Initialize position error to a large value
      double pos_err_prev = 1e6;                // Previous position error for convergence check
      Vec3d dir_prev = dir;                     // Previous direction vector for convergence check
      VecXd zf_prev = VecXd::Zero(STATE_SIZE);  // Previous state vector for convergence check

      MatXd store_mat_prev;  // Previous store matrix for convergence check

      int n_iter = 20;  // Number of iterations for coarse correction

      if (config.straight_ray == false) {
        // if straight ray, we can skip coarse correction
        for (int iter = 0; iter < n_iter; ++iter) {
          // Propagate the ray
          zf = propagate_ray(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, store_mat, false);

          pos_err = (zf.segment<3>(0) - x2).norm();  // Position error in km

          // Break the loop if coarse update no longer improves the position error
          if ((abs(pos_err - pos_err_prev) < 1e-3) || (pos_err > pos_err_prev)
              || (pos_err <= config.corr_tol / 1000)) {
            if (debug_corr) {
              std::cout << "[trace ray] Exit Coarse Update - Err: " << pos_err * 1000 << " m"
                        << std::endl;
            }

            // If the position error is larger than the previous one, restore the
            // previous direction vector
            if (pos_err > pos_err_prev) {
              if (debug_corr) {
                std::cout << "[trace ray] Position error increased, restoring "
                             "previous direction vector - Err: "
                          << pos_err_prev * 1000 << " m" << std::endl;
              }
              dir = dir_prev;              // Restore the previous direction vector
              zf = zf_prev;                // Restore the previous state vector
              store_mat = store_mat_prev;  // Restore the previous store matrix
            }

            break;  // Stop the coarse iteration if the position error is small
                    // enough
          } else {
            // The error is still large, continue with the next iteration

            // First correct the propagation time based on the TEC delay
            correct_proptime(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf);

            // Print debug info
            if (debug_corr) {
              double freq2 = freq_Hz * freq_Hz;  // Frequency squared
              double freq3 = freq2 * freq_Hz;    // Frequency cubed
              double freq4 = freq3 * freq_Hz;    // Frequency to the fourth power
              Vec3d final_pos = zf.segment<3>(0);
              double delay_tec_m = p_coeff * zf(7) / freq2;  // TEC delay in meters
              double delay_second_m = zf(8) / (freq3);
              double third_delay_m = zf(9) / (freq4);
              double C_ms = C * 1000;                     // Speed of light in m/s
              double delay_tec_sec = delay_tec_m / C_ms;  // TEC delay in seconds
              double computed_prop_time
                  = zf(6) + delay_tec_sec;                     // Final time in seconds since J2000
              double dist_prop = zf(6) * C;                    // Distance traveled
              double dist_straight = (final_pos - x1).norm();  // Straight distance
              double delay_bend_m = (dist_prop - dist_straight) * 1000;  // Bend delay in meters

              std::cout << "  Iteration: " << iter << " -----------------------------" << std::endl;
              std::cout << "    Position error  [m]: " << pos_err * 1000 << std::endl;
              std::cout << "    Time error      [s]: " << computed_prop_time - prop_time << " s  ( "
                        << C * 1000 * (computed_prop_time - prop_time) << " [m])" << std::endl;
              std::cout << "    Delay TEC       [m]: " << delay_tec_m << std::endl;
              std::cout << "    Delay bend      [m]: " << delay_bend_m << std::endl;
              std::cout << "    Delay 2nd       [m]: " << delay_second_m << std::endl;
              std::cout << "    Delay 3rd       [m]: " << third_delay_m << std::endl;
              std::cout << " " << std::endl;
            }

            // Store the previous state for convergence check
            pos_err_prev = pos_err;      // Store the previous position error
            dir_prev = dir;              // Store the previous direction vector
            zf_prev = zf;                // Store the previous state vector
            store_mat_prev = store_mat;  // Store the previous store matrix

            // Next correct the initial angle of the ray
            correct_dir_neldermead(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf, store_mat,
                                   coarse_mode, debug_corr);
          }  // End of coarse correction loop
        }
      }

      // Next Run optimization to do finer correction -----------------------
      if ((pos_err * 1000 > config.corr_tol) && (config.fine_correction)) {
        if (debug_corr) {
          std::cout << "[trace ray] Start Finer Update " << std::endl;
        }

        MatXd storemat_dummy;  // Dummy matrix for dzds
        if (correction_method == "neldermead") {
          coarse_mode = false;  // Disable coarse mode

          correct_dir_neldermead(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf,
                                 storemat_dummy, coarse_mode, debug_corr);

        } else if (correction_method == "newton") {
          correct_dir_newton(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf, debug_corr);
        } else {
          throw std::runtime_error("Unknown correction method: " + correction_method);
        }

        // Apply the final propagation time correction
        correct_proptime(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf);
      }  // Fine correction

    }  // endif apply corrections

    // compute the path profile
    if (debug_corr) {
      if (correction) {
        std::cout << "[trace ray] correction applied" << std::endl;
      } else {
        std::cout << "[trace ray] no correction applied" << std::endl;
      }
    }

    PathProfile path_profile;
    if (correction) {
      // If correction is applied, compute the path profile with the stored n and
      // grad_n
      path_profile = compute_path_profile(prop_time, dir, sf, x1, x2, epoch_utc_rx, config,
                                          store_mat, true, debug_prop);
    } else {
      // If no correction is applied, compute the path profile directly
      path_profile = compute_path_profile(prop_time, dir, sf, x1, x2, epoch_utc_rx, config,
                                          store_mat, false, debug_prop);
    }

    return path_profile;  // Return the computed path profile
  }

  void correct_proptime(double& prop_time, Vec3d& dir, double sf, const Vec3d& x1, const Vec3d& x2,
                        double t_rx, RayTraceConfig config, VecXd& zf) {
    // Unpack the configuration
    double freq_Hz = config.freq_Hz;  // Frequency in Hz

    // Update the transmission time based on the tec delay
    double t_tx_orig = t_rx - prop_time;  // Store the original transmission time
    double prop_time_orig = prop_time;    // Store the original propagation time
    Vec3d dir_orig = dir;                 // Store the original direction vector
    Vec3d tx_pos_orig = x1;

    // constants
    double C_ms = C * 1000;            // Speed of light in m/s
    double freq2 = freq_Hz * freq_Hz;  // Frequency squared
    double freq3 = freq2 * freq_Hz;    // Frequency cubed
    double freq4 = freq3 * freq_Hz;    // Frequency to the fourth power

    // Compute the final state vector
    Vec3d final_pos = zf.segment<3>(0);
    Vec3d final_dir = zf.segment<3>(3).normalized();
    double tec_delay_m = p_coeff * zf(7) / freq2;      // TEC delay in meters
    double tec_delay_sec = tec_delay_m / C_ms;         // TEC delay in s
    double second_delay_sec = zf(8) / (freq3) / C_ms;  // Second order delay in m
    double third_delay_sec = zf(9) / (freq4) / C_ms;   // Third order delay in m
    double computed_prop_time
        = zf(6) + tec_delay_sec + second_delay_sec + third_delay_sec;  // Final time in seconds

    // Update the transmission time so that prop_time = computed_prop_time
    prop_time = computed_prop_time;
  }

  void correct_dir_neldermead(double& prop_time, Vec3d& dir, double sf, const Vec3d& x1,
                              const Vec3d& x2, double t_rx, RayTraceConfig config, VecXd& zf,
                              MatXd& store_mat, bool coarse_mode, bool debug) {
    Vec3d dir_orig = dir;
    Vec2d azel_init = unitvec_to_azel(dir);

    bool use_precomputed_dzds = false;
    if (coarse_mode) use_precomputed_dzds = true;

    std::function<double(Vec2d)> func = [&](const Vec2d& azel) {
      // Convert azimuth and elevation to direction vector
      Vec3d dir_prop = azel_to_unitvec(azel(0), azel(1));
      // Perturbed ray propagation
      VecXd zf_prop = propagate_ray(prop_time, dir_prop, sf, x1, x2, t_rx, config, store_mat,
                                    use_precomputed_dzds);
      // Compute the final position and time for the perturbed ray
      Vec3d pos_prop = zf_prop.segment<3>(0);

      // Compute the position error
      double pos_err = (pos_prop - x2).norm();

      return pos_err;  // Return the position error as the objective function
    };

    // Set up the Nelder-Mead optimizer
    double step_azel, tol, fval_tol;
    bool debug_nm;
    int nm_iter;

    NelderMead nm(func);

    if (coarse_mode) {
      // coarse correction mode (reuse the precomputed dzds to save time)
      step_azel = 1e-4;  // Larger step size for coarse correction
      tol = 1e-8;        // Tolerance for convergence in coarse mode
      debug_nm = false;
      nm_iter = 200;    // Number of iterations for coarse mode
      fval_tol = 1e-6;  // terminate if the function value is below this threshold
      nm.set_initial_simplex(azel_init, step_azel,
                             debug_nm);  // Set initial simplex for coarse mode
    } else {
      // Fine correction mode (compute the refractive index at each step)
      debug_nm = debug;                            // Use the debug flag for Nelder-Mead
      double f0 = (zf.segment<3>(0) - x2).norm();  // Position error for the initial val
      if (f0 <= 0.01)                              // If the initial position error is small
        step_azel = 1e-8;                          // Smaller step size for fine correction
      else if (f0 <= 0.1)
        step_azel = 1e-7;  // Larger step size for fine correction
      else
        step_azel = 1e-6;  // Default step size for fine correction

      tol = 1e-4;                         // Tolerance for convergence in fine mode
      fval_tol = config.corr_tol / 1000;  // terminate if the function value is
                                          // below this threshold (converted to km)
      nm_iter = 50;                       // Number of iterations for fine mode
      nm.set_initial_simplex(azel_init, f0, step_azel,
                             debug_nm);  // Set initial simplex for fine mode
    }

    // Run the Nelder-Mead optimization
    Vec2d min_azel = nm.minimize(nm_iter, tol, fval_tol, debug_nm);

    // Update the propagation time based on the new azimuth and elevation
    Vec2d azel_diff = min_azel - azel_init;
    dir = azel_to_unitvec(min_azel(0), min_azel(1));

    if (debug) {
      std::cout << " " << std::endl;
      std::cout << "    [Correct Nelder Mead] " << std::endl;
      std::cout << "       Derivative precomputed: " << use_precomputed_dzds << std::endl;
      std::cout << "       Iteration used: " << nm.get_iter() << std::endl;
      std::cout << "       Min Position error [m]: " << nm.get_min_val() * 1000.0 << std::endl;
      std::cout << "       Azimuth correction [rad]: " << azel_diff(0) << std::endl;
      std::cout << "       Elevation correction [rad]: " << azel_diff(1) << std::endl;
      std::cout << "       Dir correction []:  " << (dir - dir_orig).transpose() << std::endl;
      std::cout << " " << std::endl;
    }
  }

  void correct_dir_newton(double& prop_time, Vec3d& dir, double sf, const Vec3d& x1,
                          const Vec3d& x2, double t_rx, RayTraceConfig config, VecXd& zf,
                          bool debug) {
    // Use Euler for integration to save time
    std::string integ_method = "euler";

    // Compute errors -----------------------------------------------
    Vec3d pos_err = x2 - zf.segment<3>(0);

    // Storage
    double prop_time_orig = prop_time;  // Store the original transmission time
    Vec3d dir_orig = dir;               // Store the original direction vector

    // Compute Jacobian for the inputs ---------------------------------
    Vec2d azel_orig = unitvec_to_azel(dir_orig);
    VecXd zf_pert_plus(8), zf_pert_minus(8);
    MatXd J(3, 2);        // solve for azimuth and elevation corrections
    double delta = 1e-7;  // perturbation value
    int num_cols = 2;     // Number of columns in the Jacobian (azel)

    for (int j = 0; j < num_cols; ++j) {
      // perturb the input
      Vec2d azel_pert_plus = azel_orig;   // perturbation in azimuth and elevation
      Vec2d azel_pert_minus = azel_orig;  // perturbation in azimuth and elevation

      // perturb the direction vector
      azel_pert_plus[j] = azel_orig[j] + delta;
      azel_pert_minus[j] = azel_orig[j] - delta;

      Vec3d dir_pert_plus = azel_to_unitvec(azel_pert_plus(0), azel_pert_plus(1));
      Vec3d dir_pert_minus = azel_to_unitvec(azel_pert_minus(0), azel_pert_minus(1));

      // Perturbed ray propagation (use precomputed dzds to save time)
      MatXd dzds;  // Dummy matrix for dzds
      zf_pert_plus = propagate_ray(prop_time, dir_pert_plus, sf, x1, x2, t_rx, config, dzds, false);
      zf_pert_minus
          = propagate_ray(prop_time, dir_pert_minus, sf, x1, x2, t_rx, config, dzds, false);

      // Compute the final position and time for the perturbed ray
      Vec3d pos_pert_plus = zf_pert_plus.segment<3>(0);
      Vec3d pos_pert_minus = zf_pert_minus.segment<3>(0);

      // Fill the Jacobian
      Vec3d Jcol = (pos_pert_plus - pos_pert_minus) / (2 * delta);

      if (debug) {
        std::cout << " j: " << j << std::endl;
        std::cout << "   Perturbed position (plus): " << pos_pert_plus.transpose() << std::endl;
        std::cout << "   Perturbed position (minus): " << pos_pert_minus.transpose() << std::endl;
        std::cout << "   Position error: " << (pos_pert_plus - pos_pert_minus).transpose()
                  << std::endl;
        std::cout << "   Jcol[" << j << "] = " << Jcol.transpose() << std::endl;
      }

      J.col(j) = Jcol;
    }

    // Solve the linear system J * delta = err
    Vec2d azel_corr = J.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(pos_err);

    double alpha = 1.0;  // learning rate for the correction
    Vec2d new_azel = azel_orig + alpha * azel_corr;

    // Apply the correction to the ray
    dir = azel_to_unitvec(new_azel(0), new_azel(1));

    // Correction information
    if (debug) {
      std::cout << " " << std::endl;
      std::cout << " Position error: " << pos_err.transpose() << std::endl;
      std::cout << " Jacoban matrix [3x2]: " << std::endl;
      for (int i = 0; i < J.rows(); ++i) {
        std::cout << "  " << J.row(i) << std::endl;
      }
      std::cout << " Azel correction [rad]: " << alpha * azel_corr.transpose() << std::endl;
      std::cout << " dir correction []:  " << (dir - dir_orig).transpose() << std::endl;
      std::cout << " " << std::endl;
    }
  }  // end of correct_newton

}  // namespace pecsim
