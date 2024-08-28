/**
 * @file RadioMeasurement.cpp
 * @author Stanford NAV LAB
 * @brief Class for Radionavigation measurements
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/radio_measurement.h"

#include "lupnt/measurements/comm_utils.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {
  Real ComputeOneWayRange(VecX r_tx, VecX r_rx, Real offset) {
    Real rho_rx = (r_tx - r_rx).norm();
    return rho_rx + offset;
  };

  Real ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx, Real dt_rx, Real offset) {
    // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
    Real rho_rx = (r_tx - r_rx).norm();
    Real P_rx = rho_rx + C * (dt_rx - dt_tx) + offset;
    return P_rx;
  };

  Real ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                              Real dt_rx_dot, Real offset) {
    VecX e_rx = (r_tx - r_rx).normalized();
    Real prr = e_rx.dot(v_tx - v_rx) + C * (dt_rx_dot - dt_tx_dot) + offset;
    return prr;
  };

  Real ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                           Real dt_rx_dot, Real f, Real offset) {
    Real f_D
        = -f / C * ComputePseudorangerate(r_tx, r_rx, v_tx, v_rx, dt_tx_dot, dt_rx_dot, offset);
    return f_D;
  };

  Real ComputeOneWayRangeLTR(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Real dt_tx, Real dt_rx,
                             BodyData tx_center_body, BodyData rx_center_body, bool is_bodyfixed_tx,
                             bool is_bodyfixed_rx, Real hardware_delay) {
    // Get the GM of the central body
    Real mu_tx = tx_center_body.GM;
    Real mu_rx = rx_center_body.GM;

    // transmitter and receiver states
    Vec6 xi = rv_tx;
    Vec6 x0 = rv_rx;

    // link xi->x0
    Vec3 r0 = x0.segment(0, 3);
    Vec3 v0 = x0.segment(3, 3);
    Vec3 ri = xi.segment(0, 3);
    Vec3 vi = xi.segment(3, 3);

    Vec3 a0_r0 = -mu_rx / pow(r0.norm(), 3) * r0;
    Vec3 ai_ri = -mu_tx / pow(ri.norm(), 3) * ri;

    // hardware delays
    Real tau_d_rx = hardware_delay;  // receiver delay for downlink (xi->x0)

    // solve for tau_d (downlink time)
    int max_iter = 5;
    Real tau_d = 0.0;
    Real tau_d_prev = 0.0;

    Vec3 r0_p, rid_p, rho_ad;
    Real r0_p_norm, rid_p_norm, rho_ad_norm;

    // Solve for Downlink
    for (int i = 0; i < max_iter; i++) {
      Real rx_t = epoch_rx - tau_d_rx;
      Real tx_t = epoch_rx - tau_d_rx - tau_d;

      if (!is_bodyfixed_rx) {  // is not ground station
        r0_p = r0 - v0 * tau_d_rx + 1 / 2 * a0_r0 * pow(tau_d_rx, 2);
      } else {
        // Fixed to surface frame -> convert to ITRF frame
        r0_p = ConvertFrame(rx_t, r0, rx_center_body.fixed_frame, Frame::ITRF);
      }
      if (!is_bodyfixed_tx) {  // is not ground station
        rid_p = ri - vi * (tau_d_rx + tau_d) + 1 / 2 * ai_ri * pow(tau_d_rx + tau_d, 2);
      } else {
        // Fixed to surface frame -> convert to ITRF frame
        rid_p = ConvertFrame(tx_t, ri, tx_center_body.fixed_frame, Frame::ITRF);
      }
      rho_ad = rid_p - r0_p;

      // norms
      r0_p_norm = r0_p.norm();
      rid_p_norm = rid_p.norm();
      rho_ad_norm = rho_ad.norm();

      // Compoensate for relativistic effects (Shapiro time delay)
      // double shapiro = 2 * mu_rx / C *
      //                  log((r0_p_norm + rid_p_norm + rho_ad_norm) /
      //                      (r0_p_norm + rid_p_norm - rho_ad_norm));

      tau_d = rho_ad.norm();  // + shapiro;

      if (fabs(tau_d.val() - tau_d_prev.val()) < 1e-13) {  // goes under pico-second
        break;
      } else {
        tau_d_prev = tau_d;
      }
    }

    Real rho_d = C * tau_d + (dt_rx - dt_tx) * C;

    return rho_d;
  };

  Real ComputeTwoWayRangeLTR(Real epoch_rx, Vec6 rv_target_tr, Vec6 rv_rx_tr,
                             BodyData target_center_body, BodyData rx_center_body,
                             bool is_bodyfixed_target, bool is_bodyfixed_rx, Real hardware_delay,
                             Real additional_delay) {
    // hardware delays
    Real tau_d_rx = hardware_delay;  // receiver delay for downlink (xi->x0)
    Real tau_u_rx = hardware_delay;  // receiver delay for uplink (x0->xi)
    Real tau_d_tx = hardware_delay;  // transmitter delay for downlink (xi->x0)

    // solve for tau_d (downlink time, target->rx)
    Real rho_d = ComputeOneWayRangeLTR(epoch_rx, rv_target_tr, rv_rx_tr, 0.0, 0.0,
                                       target_center_body, rx_center_body, is_bodyfixed_target,
                                       is_bodyfixed_rx, tau_d_rx + additional_delay);
    Real tau_d = rho_d / C;

    // solve for tau_u (uplink time, rx->target)
    Real tau_c_pp = tau_u_rx + tau_d_tx + tau_d_rx;  // total hardware delay
    Real delay_uplink = tau_c_pp + tau_d;            // total delay for uplink w.r.t epoch_rx
    Real rho_u = ComputeOneWayRangeLTR(epoch_rx, rv_rx_tr, rv_target_tr, 0.0, 0.0, rx_center_body,
                                       target_center_body, is_bodyfixed_rx, is_bodyfixed_target,
                                       delay_uplink + additional_delay);
    Real tau_u = rho_u / C;

    Real rho_ud = C / 2 * (tau_u + tau_d);

    return rho_ud;
  };

  Real ComputeOneWayRangeRateLTR(Real epoch_rx, Vec6 rv_tx_tr, Vec6 rv_rx_tr, Real dt_dot_tx,
                                 Real dt_dot_rx, BodyData target_center_body,
                                 BodyData rx_center_body, bool is_bodyfixed_target,
                                 bool is_bodyfixed_rx, Real hardware_delay, double T_I) {
    Real rho_d = ComputeOneWayRangeLTR(epoch_rx, rv_tx_tr, rv_rx_tr, 0.0, 0.0, target_center_body,
                                       rx_center_body, is_bodyfixed_target, is_bodyfixed_rx,
                                       hardware_delay);
    Real rho_d_past = ComputeOneWayRangeLTR(epoch_rx, rv_tx_tr, rv_rx_tr, 0.0, 0.0,
                                            target_center_body, rx_center_body, is_bodyfixed_target,
                                            is_bodyfixed_rx, hardware_delay + T_I);

    Real rho_dot = (rho_d - rho_d_past) / T_I + C * (dt_dot_rx - dt_dot_tx);

    return rho_dot;
  }

  Real ComputeTwoWayRangeRateLTR(Real epoch_rx, Vec6 rv_target_tr, Vec6 rv_rx_tr,
                                 BodyData target_center_body, BodyData rx_center_body,
                                 bool is_bodyfixed_target, bool is_bodyfixed_rx,
                                 Real hardware_delay, double T_I) {
    Real rho_ud = ComputeTwoWayRangeLTR(epoch_rx, rv_target_tr, rv_rx_tr, target_center_body,
                                        rx_center_body, is_bodyfixed_target, is_bodyfixed_rx,
                                        hardware_delay, 0);
    Real rho_ud_past = ComputeTwoWayRangeLTR(epoch_rx, rv_target_tr, rv_rx_tr, target_center_body,
                                             rx_center_body, is_bodyfixed_target, is_bodyfixed_rx,
                                             hardware_delay, T_I);

    Real rho_dot = (rho_ud - rho_ud_past) / T_I;

    return rho_dot;
  }

  double ComputePnRangeErrorCTL(double PRC_N0, double B_L, double Tc, Modulation modulation_type) {
    double sigma = 0.0;
    double f_RC = 1 / (2 * Tc);

    // Thermal noise
    sigma = 1 / sqrt(2) * C / (8 * f_RC) * sqrt(B_L / PRC_N0);

    // Error degrade for GMSK + PN

    return sigma;
  }

  double ComputePnRangeErrorOL(double PRC_N0, double TI, double Tc, Modulation modulation_type) {
    double sigma = 0.0;
    double f_RC = 1 / (2 * Tc);

    // Thermal noise
    sigma = 1 / sqrt(32 * PI * PI) * (C / f_RC) * sqrt(1 / PRC_N0 / TI);

    // Error degrade for GMSK + PN

    return sigma;
  }

  double ComputeRangeRateErrorOneWay(double B_L_carrier, double f_C, double T_s, double T_I,
                                     double PT_N0, double sigma_y_1s, Modulation modulation_type,
                                     double m_R) {
    // Thermal noise
    double rho_L = ComputeCarrierLoopSNR(PT_N0, B_L_carrier, T_s, modulation_type, m_R);
    double sigma_vn = sqrt(2 / rho_L) * C / (2 * PI * f_C * T_I);

    // phase noise contribution
    double sigma_y_T = sigma_y_1s / sqrt(T_I);
    double sigma_vf = C * sigma_y_T;

    // phase scintillation
    double sigma_vs = 0.0;  // Asssume 0

    // Total Doppler Error
    double sigma_v = sqrt(pow(sigma_vn, 2) + pow(sigma_vf, 2) + pow(sigma_vs, 2));

    return sigma_v;
  }

  double ComputeRangeRateErrorTwoWay(double B_L_carrier, double f_C, double T_s, double T_I,
                                     double PT_N0, double sigma_y_1s, double G,
                                     Modulation modulation_type, double m_R) {
    // Thermal noise
    double rho_L = ComputeCarrierLoopSNR(PT_N0, B_L_carrier, T_s, modulation_type, m_R);
    double sigma_vnu = sqrt(1 / 2) * (C / (2 * PI * f_C * T_I)) * G / sqrt(rho_L);
    double sigma_vnd = sqrt(2 / rho_L) * C / (2 * PI * f_C * T_I) / sqrt(rho_L);

    double sigma_vn = sqrt(pow(sigma_vnu, 2) + pow(sigma_vnd, 2));

    // phase noise contribution
    double sigma_y_T = sigma_y_1s / sqrt(T_I);
    double sigma_vf = C * sigma_y_T / sqrt(2);

    // phase scintillation
    double sigma_vs = 0.0;  // Asssume 0

    // Total Doppler Error
    double sigma_v = sqrt(pow(sigma_vn, 2) + pow(sigma_vf, 2) + pow(sigma_vs, 2));

    // Error degrade for GMSK + PN

    return sigma_v;
  }

}  // namespace lupnt
