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

namespace lupnt {
Real RadioMeasurement::ComputeOneWayRange(VecX r_tx, VecX r_rx, Real offset) {
  Real rho_rx = (r_tx - r_rx).norm();
  return rho_rx + offset;
};

Real RadioMeasurement::ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx,
                                          Real dt_rx, Real offset) {
  // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
  Real rho_rx = (r_tx - r_rx).norm();
  Real P_rx = rho_rx + C * (dt_rx - dt_tx) + offset;
  return P_rx;
};

Real RadioMeasurement::ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx,
                                              VecX v_rx, Real dt_tx_dot,
                                              Real dt_rx_dot, Real offset) {
  VecX e_rx = (r_tx - r_rx).normalized();
  Real prr = e_rx.dot(v_tx - v_rx) + C * (dt_rx_dot - dt_tx_dot) + offset;
  return prr;
};

Real RadioMeasurement::ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx,
                                           VecX v_rx, Real dt_tx_dot,
                                           Real dt_rx_dot, Real f,
                                           Real offset) {
  Real f_D = -f / C *
             ComputePseudorangerate(r_tx, r_rx, v_tx, v_rx, dt_tx_dot,
                                    dt_rx_dot, offset);
  return f_D;
};

Real RadioMeasurement::ComputeOneWayRangeLTR(VecX rv_tx, VecX rv_rx, Real dt_tx,
                                             Real dt_rx, double mu,
                                             Real hardware_delay) {
  // transmitter and receiver states
  VecX xi = rv_tx;
  VecX x0 = rv_rx;

  // link x0->xi
  VecX r0 = x0.segment(0, 3);
  VecX v0 = x0.segment(3, 3);
  VecX ri = xi.segment(0, 3);
  VecX vi = xi.segment(3, 3);

  VecX a0_r0 = -mu / pow(r0.norm(), 3) * r0;
  VecX ai_ri = -mu / pow(ri.norm(), 3) * ri;

  // hardware delays
  Real tau_d_rx = hardware_delay;  // receiver delay for downlink (xi->x0)

  // solve for tau_d (downlink time)
  int max_iter = 5;
  Real tau_d = 0.0;
  Real tau_d_prev = 0.0;

  VecX r0_p(3), rid_p(3), rho_ad(3);
  Real r0_p_norm, rid_p_norm, rho_ad_norm;

  for (int i = 0; i < max_iter; i++) {
    r0_p = r0 - v0 * tau_d_rx + 1 / 2 * a0_r0 * pow(tau_d_rx, 2);
    rid_p =
        ri - vi * (tau_d_rx + tau_d) + 1 / 2 * ai_ri * pow(tau_d_rx + tau_d, 2);
    rho_ad = rid_p - r0_p;

    // norms
    r0_p_norm = r0_p.norm();
    rid_p_norm = rid_p.norm();
    rho_ad_norm = rho_ad.norm();

    // Compoensate for relativistic effects (Shapiro time delay)
    tau_d = rho_ad.norm() + 2 * mu / C *
                                log((r0_p_norm + rid_p_norm + rho_ad_norm) /
                                    (r0_p_norm + rid_p_norm - rho_ad_norm));
    if (fabs(tau_d.val() - tau_d_prev.val()) <
        1e-13) {  // goes under pico-second
      break;
    } else {
      tau_d_prev = tau_d;
    }
  }

  Real rho_d = C * tau_d + (dt_rx - dt_tx) * C;

  return rho_d;
};

Real RadioMeasurement::ComputeTwoWayRangeLTR(VecX rv_target_tr, VecX rv_rx_tr,
                                             double mu, Real hardware_delay) {
  // transmitter and receiver states
  VecX xi = rv_target_tr;
  VecX x0 = rv_rx_tr;

  // link x0->xi
  VecX r0 = x0.segment(0, 3);
  VecX v0 = x0.segment(3, 3);
  VecX ri = xi.segment(0, 3);
  VecX vi = xi.segment(3, 3);

  VecX a0_r0 = -mu / pow(r0.norm(), 3) * r0;
  VecX ai_ri = -mu / pow(ri.norm(), 3) * ri;

  // hardware delays
  Real tau_d_rx = hardware_delay;  // receiver delay for downlink (xi->x0)
  Real tau_u_rx = hardware_delay;  // receiver delay for uplink (x0->xi)
  Real tau_d_tx = hardware_delay;  // transmitter delay for downlink (xi->x0)

  // solve for tau_d (downlink time) -----------------------------------
  int max_iter = 5;
  Real tau_d = 0.0;
  Real tau_d_prev = 0.0;

  VecX r0_p(3), rid_p(3), rho_ad(3);
  Real r0_p_norm, rid_p_norm, rho_ad_norm;

  for (int i = 0; i < max_iter; i++) {
    r0_p = r0 - v0 * tau_d_rx + 1 / 2 * a0_r0 * pow(tau_d_rx, 2);
    rid_p =
        ri - vi * (tau_d_rx + tau_d) + 1 / 2 * ai_ri * pow(tau_d_rx + tau_d, 2);
    rho_ad = rid_p - r0_p;

    // norms
    r0_p_norm = r0_p.norm();
    rid_p_norm = rid_p.norm();
    rho_ad_norm = rho_ad.norm();

    // Compoensate for relativistic effects (Shapiro time delay)
    tau_d = rho_ad.norm() + 2 * mu / C *
                                log((r0_p_norm + rid_p_norm + rho_ad_norm) /
                                    (r0_p_norm + rid_p_norm - rho_ad_norm));
    if (fabs(tau_d.val() - tau_d_prev.val()) <
        1e-13) {  // goes under pico-second
      break;
    } else {
      tau_d_prev = tau_d;
    }
  }

  // solve for tau_u (uplink time) -----------------------------------
  Real tau_u = 0.0;
  Real tau_u_prev = 0.0;
  Real tau_c_pp = tau_u_rx + tau_d_tx + tau_d_rx;
  max_iter = 5;
  VecX r0u_pp(3), ri_pp(3), rho_au(3);
  Real r0u_pp_norm, ri_pp_norm, rho_au_norm;

  for (int i = 0; i < max_iter; i++) {
    ri_pp =
        ri - vi * (tau_c_pp + tau_d) + 1 / 2 * ai_ri * pow(tau_c_pp + tau_d, 2);
    r0u_pp = r0 - v0 * (tau_c_pp + tau_d + tau_u) +
             1 / 2 * a0_r0 * pow(tau_c_pp + tau_d + tau_u, 2);

    rho_au = ri_pp - r0u_pp;

    // norms
    ri_pp_norm = ri_pp.norm();
    r0u_pp_norm = r0u_pp.norm();
    rho_au_norm = rho_au.norm();

    // Compoensate for relativistic effects (Shapiro time delay)
    tau_u = rho_au.norm() + 2 * mu / C *
                                log((ri_pp_norm + r0u_pp_norm + rho_au_norm) /
                                    (ri_pp_norm + r0u_pp_norm - rho_au_norm));

    if (fabs(tau_u.val() - tau_u_prev.val()) <
        1e-13) {  // goes under pico-second
      break;
    } else {
      tau_u_prev = tau_u;
    }
  }

  Real rho_ud = C / 2 * (tau_u + tau_d);

  return rho_ud;
};

double RadioMeasurement::ComputePnRangeErrorCTL(double PRC_N0, double B_L,
                                                double Tc) {
  double sigma = 0.0;
  double f_RC = 1 / (2 * Tc);

  sigma = 1 / sqrt(2) * C / (8 * f_RC) * sqrt(B_L / PRC_N0);

  return sigma;
}

double RadioMeasurement::ComputePnRangeErrorOL(double PRC_N0, double TI,
                                               double Tc) {
  double sigma = 0.0;
  double f_RC = 1 / (2 * Tc);

  sigma = 1 / sqrt(32 * PI * PI) * (C / f_RC) * sqrt(1 / PRC_N0 / TI);

  return sigma;
}

double ComputeRangeRateErrorOneWay(double B_L_carrier, double f_C, double T_s,
                                   double T_I, double PT_N0, double sigma_y_1s,
                                   Modulation carrier_type = Modulation::BPSK) {
  　  // Thermal noise
      double rho_L =
          ComputeCarrierLoopSNR(PT_N0, B_L_carrier, T_s, carrier_type);
  double sigma_vn = sqrt(2 / rho_L) * C / (2 * PI * f_C * T_I);

  // phase noise contribution
  double sigma_y_T = sigma_y_1s / sqrt(T_s);
  double sigma_vf = C * sigma_y_T;

  // phase scintillation
  double sigma_vs = 0.0;  // Asssume 0

  // Total Doppler Error
  double sigma_v = sqrt(pow(sigma_vn, 2) + pow(sigma_vf, 2) + pow(sigma_vs, 2));
  return sigma_v;
}

double ComputeRangeRateErrorTwoWay(double B_L_carrier, double f_C, double T_s,
                                   double T_I, double PT_N0, double sigma_y_1s,
                                   double G,
                                   Modulation carrier_type = Modulation::BPSK) {
  　  // Thermal noise
      double rho_L =
          ComputeCarrierLoopSNR(PT_N0, B_L_carrier, T_s, carrier_type);
  double sigma_vnu =
      sqrt(1 / 2) * (C / (2 * PI * f_C * T_I)) * pow(G) / sqrt(rho_L);
  double sigma_vnd = sqrt(2 / rho_L) * C / (2 * PI * f_C * T_I) / sqrt(rho_L);

  double sigma_vn = sqrt(pow(sigma_vnu, 2) + pow(sigma_vnd, 2));

  // phase noise contribution
  double sigma_y_T = sigma_y_1s / sqrt(T_s);
  double sigma_vf = C * sigma_y_T / sqrt(2);

  // phase scintillation
  double sigma_vs = 0.0;  // Asssume 0

  // Total Doppler Error
  double sigma_v = sqrt(pow(sigma_vn, 2) + pow(sigma_vf, 2) + pow(sigma_vs, 2));
  return sigma_v;
}

double GetTransponderTurnAroundRatio(FrequencyBand fbu, FrequencyBand fbd) {
  // https://deepspace.jpl.nasa.gov/dsndocs/810-005/201/201B.pdf
  double G = 1.0;
  // fbu: S, X, Ka,  fbd: fbu: S, X, Ka
  if (fbu == FrequencyBand::S && fbd == FrequencyBand::S) {
    G = 240 / 221;
  } else if (fbu == FrequencyBand::S && fbd == FrequencyBand::X) {
    G = 880 / 221;
  } else if (fbu == FrequencyBand::S && fbd == FrequencyBand::Ka) {
    G = 15.071;
  } else if (fbu == FrequencyBand::X && fbd == FrequencyBand::S) {
    G = 240 / 749;
  } else if (fbu == FrequencyBand::X && fbd == FrequencyBand::X) {
    G = 880 / 749;
  } else if (fbu == FrequencyBand::X && fbd == FrequencyBand::Ka) {
    G = 4.4506;
  } else if (fbu == FrequencyBand::Ka && fbd == FrequencyBand::S) {
    G = 0.066959;
  } else if (fbu == FrequencyBand::Ka && fbd == FrequencyBand::X) {
    G = 0.24561;
  } else if (fbu == FrequencyBand::Ka && fbd == FrequencyBand::Ka) {
    G = 0.92982;
  }

  return G;
}

double ComputeCarrierLoopSNR(double PT_N0, double B_L_carrier, double T_s,
                             Modulation carrier_type) {
  double S_L = 1.0;
  double EsN0 = PT_N0 * T_s;

  // Compute carier loop signal-to-noise ratio
  switch (carrier_type) {
    case Modulation::BPSK:
      S_L = 2 * EsN0 / (1 + 2 * EsN0);
      break;
    case Modulation::QPSK:
      double tmp1 = (9 / 4) / EsN0;
      double tmp2 = (3 / 2) / pow(EsN0, 2);
      double tmp3 = (3 / 16) / pow(EsN0, st3);

      S_L = 1 / (1 + tmp1 + tmp2 + tmp3);
      break;
    case Modulation::OQPSK:
      double tmp1 = (9 / 4) / EsN0;
      double tmp2 = (3 / 2) / pow(EsN0, 2);
      double tmp3 = (3 / 16) / pow(EsN0, 3);

      S_L = 1 / 4 / (1 + tmp1 + tmp2 + tmp3);
      break;
    case Modulation::GMSK:
      // same as OQPSK
      // Reference: GMSK Modulation for Deep Space Application
      // https://ieeexplore-ieee-org.stanford.idm.oclc.org/stamp/stamp.jsp?tp=&arnumber=6187097
      double tmp1 = (9 / 4) / EsN0;
      double tmp2 = (3 / 2) / pow(EsN0, 2);
      double tmp3 = (3 / 16) / pow(EsN0, 3);

      S_L = 1 / 4 / (1 + tmp1 + tmp2 + tmp3);
      break;
  }
  double rho_L = PT_N0 * S_L / B_L_carrier;

  return rho_L;
}

}  // namespace lupnt