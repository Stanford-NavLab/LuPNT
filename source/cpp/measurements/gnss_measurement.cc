/**
 * @file gnss_measurement.cpp
 * @author Stanford NAV LAB
 * @brief Class that constructs GPS measurement data from channel observables
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "lupnt/measurements/gnss_measurement.h"

#include "lupnt/measurements/radio_measurement.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  GnssMeasurement::GnssMeasurement(const std::vector<GnssTransmission> trans)
      : trans_store(trans),
        dt_tx(trans.size()),
        dt_tx_dot(trans.size()),
        I_rx(trans.size()),
        T_rx(trans.size()),
        N_rx(trans.size()),
        eps_P(trans.size()),
        f(trans.size()),
        r_tx(3, trans.size()),
        v_tx(3, trans.size()),
        r_rx(3),
        v_rx(3),
        CN0(trans.size()),
        vis_earth(trans.size()),
        vis_moon(trans.size()),
        vis_antenna(trans.size()),
        vis_atmos(trans.size()),
        vis_ionos(trans.size()),
        rho_rx(trans.size()),
        lambda_(trans.size()),
        P_rx(trans.size()) {
    n_meas = trans.size();

    // Iterate over received transmissions and update dt_tx, I_rx, and T_rx
    int i = 0;
    for (auto tr : trans) {
      // Clock time
      t_rx = tr.t_rx;
      t_tx = tr.t_tx;

      // Clock offset
      dt_rx = tr.dt_rx;
      dt_tx[i] = tr.dt_tx;
      dt_rx_dot = tr.dt_rx_dot;
      dt_tx_dot[i] = tr.dt_tx_dot;

      // Position
      r_rx = tr.r_rx;
      r_tx.col(i) = tr.r_tx;

      // Velocity
      v_rx = tr.v_rx;
      v_tx.col(i) = tr.v_tx;

      // Propagation
      I_rx[i] = tr.I_rx;
      T_rx[i] = tr.T_rx;
      CN0[i] = tr.CN0;
      eps_P[i] = 0.0;

      // Integer Ambiguity
      N_rx[i] = 0.0;  // Todo: Add integer ambiguity models

      // Visibility
      vis_earth[i] = tr.vis_earth;
      vis_moon[i] = tr.vis_moon;
      vis_antenna[i] = tr.vis_antenna;
      vis_atmos[i] = tr.vis_atmos;
      vis_ionos[i] = tr.vis_ionos;

      f[i] = tr.freq;
      lambda_[i] = c / f[i];

      gnssr_param = tr.gnssr_param;
      chip_rate = tr.chip_rate;

      ID_tx.push_back(tr.ID_tx);

      i++;
    }
  }

  /**
   * @brief Create another measurement class instance with certain band
   *
   * @param freq_label
   * @return GnssMeasurement
   */
  GnssMeasurement GnssMeasurement::ExtractSignal(std::string freq_label) {
    std::vector<GnssTransmission> transmissions_freq;
    for (auto &tx : trans_store) {
      if (tx.freq_label == freq_label) {
        transmissions_freq.push_back(tx);
      }
    }

    return GnssMeasurement(transmissions_freq);
  }

  /***********************************************************
   * General Methods for computing Measurements
   ***********************************************************/

  VecX GnssMeasurement::ComputeGnssPseudorange(VecX r_rx, Real dt_rx, bool with_noise, int seed) {
    // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
    VecX P_rx(r_tx.cols());

    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < r_tx.cols(); i++) {
      P_rx(i) = ComputePseudorange(r_rx, r_tx.col(i), dt_tx[i], dt_rx, I_rx(i) + T_rx(i));
      if (with_noise) {
        double sigma = ComputeGnssPseudorangeNoise(CN0(i));
        P_rx(i) += sigma * distribution(generator);
      }
    }
    return P_rx;
  }

  VecX GnssMeasurement::ComputeGnssPseudorangerate(VecX r_rx, VecX v_rx, Real dt_rx_dot,
                                                   bool with_noise, int seed) {
    VecX P_rrx(r_tx.cols());

    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<double> distribution(0.0, 1.0);

    Real offset = 0;  // set pseudorange rate offset to zero

    for (int i = 0; i < r_tx.cols(); i++) {
      P_rrx(i) = ComputePseudorangerate(r_tx.col(i), r_rx, v_tx.col(i), v_rx, dt_tx_dot[i],
                                        dt_rx_dot, offset);
      if (with_noise) {
        double sigma = ComputeGnssPseudorangerateNoise(CN0(i), lambda_(i));
        P_rrx(i) += sigma * distribution(generator);
      }
    }

    return P_rrx;
  }

  VecX GnssMeasurement::ComputeGnssCarrierPhase(VecX r_rx, Real dt_rx, VecX N_rx, bool with_noise,
                                                int seed) {
    // phi_rx = c / lambda_ * (t_rx - t_tx) + c / lambda_ * (dt_rx(t_rx) -
    // dt_tx(t_tx)) + phi_rx_0 - phi_0 + N_rx + eps_phi
    VecX phi_rx(r_tx.cols());

    // check if N_rx has the same size as r_tx
    if (N_rx.size() != r_tx.cols()) {
      std::cout << "The integer ambiguity N_rx has different size than r_tx" << std::endl;
      return phi_rx;
    }

    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<double> distribution(0.0, 1.0);

    Real pr, phase;

    for (int i = 0; i < r_tx.cols(); i++) {
      pr = ComputePseudorange(r_rx, r_tx.col(i), dt_tx[i], dt_rx,
                              -I_rx(i) + T_rx(i));  // Ionosphere acts negative on phase
      phase = pr / lambda_(i) + N_rx(i);            // phase = pr/lambda + N (N=integer ambiguity)

      if (with_noise) {
        double sigma = ComputeGnssCarrierPhaseNoise(CN0(i), lambda_(i));
        phi_rx(i) += sigma * distribution(generator);
      }
    }
    return phi_rx;
  }

  /***********************************************************
   *  Methods for true measurement generation
   ***********************************************************/

  VecX GnssMeasurement::GetGnssMeasurement(std::vector<GnssMeasurementType> meas_type,
                                           bool with_noise, int seed) {
    VecX z(n_meas * meas_type.size());
    int i = 0;
    for (auto type : meas_type) {
      switch (type) {
        case GnssMeasurementType::PR:
          z.segment(i * n_meas, n_meas) = GetPseudorange(with_noise, seed);
          break;
        case GnssMeasurementType::PRR:
          z.segment(i * n_meas, n_meas) = GetPseudorangerate(with_noise, seed);
          break;
        case GnssMeasurementType::CP:
          z.segment(i * n_meas, n_meas) = GetCarrierPhase(with_noise, seed);
          break;
        default: std::cout << "Measurement type: " << type << " not supported" << std::endl; break;
      }
      i++;
    }
    return z;
  }

  VecX GnssMeasurement::GetPseudorange(bool with_noise, int seed) {
    return ComputeGnssPseudorange(r_rx, dt_rx, with_noise, seed);
  }

  VecX GnssMeasurement::GetPseudorangerate(bool with_noise, int seed) {
    return ComputeGnssPseudorangerate(r_rx, v_rx, dt_rx_dot, with_noise, seed);
  }

  VecX GnssMeasurement::GetCarrierPhase(bool with_noise, int seed) {
    return ComputeGnssCarrierPhase(r_rx, dt_rx, N_rx, with_noise, seed);
  }

  /***********************************************************
   *  Methods for predicted measurement generation
   ***********************************************************/
  VecX GnssMeasurement::GetPredictedGnssMeasurement(double epoch, Vec6 rv_pred, Vec2 clk_pred,
                                                    VecX N_pred, MatXd &H_gnss,
                                                    std::vector<GnssMeasurementType> meas_type,
                                                    Frame frame_in) {
    // number of measurements
    int n_meas_all = n_meas * meas_type.size();
    bool use_cp = false;
    for (auto type : meas_type) {
      if (type == GnssMeasurementType::CP) use_cp = true;
    }
    // determine state size
    int state_size = 8;  // position, velocity, clock
    if (use_cp) state_size = 9;

    VecX z(n_meas_all);
    H_gnss.resize(n_meas_all, state_size);
    H_gnss = MatXd::Zero(n_meas_all, state_size);
    MatXd H_pr(n_meas, state_size), H_prr(n_meas, state_size), H_cp(n_meas, state_size);

    int i = 0;

    for (auto type : meas_type) {
      switch (type) {
        case GnssMeasurementType::PR:
          H_pr = MatXd::Zero(n_meas, state_size);
          z.segment(i * n_meas, n_meas)
              = GetPredictedPseudorange(epoch, rv_pred, clk_pred, H_pr, frame_in);
          H_gnss.block(i * n_meas, 0, n_meas, 8) = H_pr;
          break;
        case GnssMeasurementType::PRR:
          H_prr = MatXd::Zero(n_meas, state_size);
          z.segment(i * n_meas, n_meas)
              = GetPredictedPseudorangerate(epoch, rv_pred, clk_pred, H_prr, frame_in);
          H_gnss.block(i * n_meas, 0, n_meas, 8) = H_prr;
          break;
        case GnssMeasurementType::CP:
          H_cp = MatXd::Zero(n_meas, state_size);
          z.segment(i * n_meas, n_meas)
              = GetPredictedCarrierPhase(epoch, rv_pred, clk_pred, N_pred, H_cp, frame_in);
          H_gnss.block(i * n_meas, 0, n_meas, 9) = H_cp;
          break;
        default: std::cout << "Measurement type: " << type << " not supported" << std::endl; break;
      }
      i++;
    }
    return z;
  }

  VecX GnssMeasurement::GetPredictedPseudorange(double epoch, Vec6 rv_pred, Vec2 clk_pred,
                                                MatXd &H_pr, Frame frame_in) {
    auto func = [epoch, frame_in, this](const Vec6 rv_in, const Vec2 clk) {
      Vec6 rv_gcrf = ConvertFrame(epoch, rv_in, frame_in, Frame::GCRF);
      Vec3 r_rx = rv_gcrf.head(3);
      Real dt_rx = clk(0);
      return ComputeGnssPseudorange(r_rx, dt_rx);
    };

    VecX z_pr_pred(n_meas);

    // break the computational graph relations before taking jacobian
    Vec6 rv_in_tmp = rv_pred.cast<double>();
    Vec2 clk_tmp = clk_pred.cast<double>();

    H_pr = jacobian(func, wrt(rv_in_tmp, clk_tmp), at(rv_in_tmp, clk_tmp), z_pr_pred);
    return z_pr_pred;
  }

  VecX GnssMeasurement::GetPredictedPseudorangeAnalyticalJacobian(double epoch, Vec6 rv_pred,
                                                                  Vec2 clk_pred, MatXd &H_pr,
                                                                  Frame frame_in) {
    auto rv_gcrf = ConvertFrame(epoch, rv_pred, frame_in, Frame::GCRF);
    Vec3 r_rx = rv_gcrf.head(3);
    Real dt_rx = clk_pred(0);

    // compute range
    VecX P_rx(r_tx.cols());
    H_pr = VecXd::Zero(r_tx.cols(), 8);

    for (int i = 0; i < r_tx.cols(); i++) {
      Real offset = I_rx(i) + T_rx(i);
      VecX r_tx_col = r_tx.col(i);
      Real dt_tx_col = dt_tx[i];
      Real rho_rx = (r_tx_col - r_rx).norm();

      // directly compute jacobian
      P_rx(i) = rho_rx + C * (dt_rx - dt_tx_col) + offset;
      H_pr(i, 0) = ((r_rx(0) - r_tx_col(0)) / rho_rx).val();
      H_pr(i, 1) = ((r_rx(1) - r_tx_col(1)) / rho_rx).val();
      H_pr(i, 2) = ((r_rx(2) - r_tx_col(2)) / rho_rx).val();
      H_pr(i, 6) = C;
    }

    return P_rx;
  }

  VecX GnssMeasurement::GetPredictedPseudorangerate(double epoch, Vec6 rv_pred, Vec2 clk_pred,
                                                    MatXd &H_prr, Frame frame_in) {
    auto func = [epoch, frame_in, this](const Vec6 rv_in, const Vec2 clk) {
      Vec6 rv_gcrf = ConvertFrame(epoch, rv_in, frame_in, Frame::GCRF);
      Vec3 r_rx = rv_gcrf.head(3);
      Vec3 v_rx = rv_gcrf.tail(3);
      Real dt_rx_dot = clk(1);
      return ComputeGnssPseudorangerate(r_rx, v_rx, dt_rx_dot);
    };

    VecX z_prr_pred(n_meas);

    // break the computational graph relations before taking jacobian
    Vec6 rv_in_tmp = rv_pred.cast<double>();
    Vec2 clk_tmp = clk_pred.cast<double>();

    H_prr = jacobian(func, wrt(rv_in_tmp, clk_tmp), at(rv_in_tmp, clk_tmp), z_prr_pred);
    return z_prr_pred;
  }

  VecX GnssMeasurement::GetPredictedCarrierPhase(double epoch, Vec6 rv_pred, Vec2 clk_pred,
                                                 VecX N_pred, MatXd &H_cp, Frame frame_in) {
    auto func = [epoch, frame_in, this](const Vec6 rv_in, const Vec2 clk, const VecX N_pred) {
      Vec6 rv_gcrf = ConvertFrame(epoch, rv_in, frame_in, Frame::GCRF);
      Vec3 r_rx = rv_gcrf.head(3);
      Real dt_rx = clk(0);
      return ComputeGnssCarrierPhase(r_rx, dt_rx, N_pred);
    };

    VecX z_cp_pred(n_meas);

    // break the computational graph relations before taking jacobian
    Vec6 rv_in_tmp = rv_pred.cast<double>();
    Vec2 clk_tmp = clk_pred.cast<double>();
    VecX N_pred_tmp = N_pred.cast<double>();

    H_cp = jacobian(func, wrt(rv_in_tmp, clk_tmp, N_pred_tmp), at(rv_in_tmp, clk_tmp, N_pred_tmp),
                    z_cp_pred);
    return z_cp_pred;
  }

  /********************
   * Noise Models
   *********************/
  VecXd GnssMeasurement::GetGnssNoiseStdVec(std::vector<GnssMeasurementType> meas_type) {
    VecXd noise(n_meas * meas_type.size());
    int i = 0;
    for (auto type : meas_type) {
      switch (type) {
        case GnssMeasurementType::PR:
          noise.segment(i * n_meas, n_meas) = GetPseudorangeNoiseStdVec();
          break;
        case GnssMeasurementType::PRR:
          noise.segment(i * n_meas, n_meas) = GetPseudorangeRateNoiseStdVec();
          break;
        case GnssMeasurementType::CP:
          noise.segment(i * n_meas, n_meas) = GetCarrierPhaseNoiseStdVec();
          break;
        default: std::cout << "Measurement type: " << type << " not supported" << std::endl; break;
      }
      i++;
    }
    return noise;
  }

  VecXd GnssMeasurement::GetPseudorangeNoiseStdVec() {
    int n_meas = r_tx.cols();
    VecXd noise(n_meas);

    for (int i = 0; i < n_meas; i++) {
      noise(i) = ComputeGnssPseudorangeNoise(CN0(i));
    }
    return noise;
  }

  VecXd GnssMeasurement::GetPseudorangeRateNoiseStdVec() {
    int n_meas = r_tx.cols();
    VecXd noise(n_meas);

    for (int i = 0; i < n_meas; i++) {
      noise(i) = ComputeGnssPseudorangerateNoise(CN0(i), lambda_(i));
    }
    return noise;
  }

  VecXd GnssMeasurement::GetCarrierPhaseNoiseStdVec() {
    int n_meas = r_tx.cols();
    VecXd noise(n_meas);

    for (int i = 0; i < n_meas; i++) {
      noise(i) = ComputeGnssCarrierPhaseNoise(CN0(i), lambda_(i));
    }
    return noise;
  }

  double GnssMeasurement::ComputeGnssPseudorangeNoise(double CN0_dB) {
    // thermal noise in DLL
    double sigma = 0.0;
    double CN0 = pow(10, CN0_dB / 10);

    // extract gnss receiver parameters
    double Bn = gnssr_param.Bn;
    double Rc = chip_rate;
    double Bfe = gnssr_param.b * Rc;
    double T = gnssr_param.T;
    double D = gnssr_param.D;
    double Tc = 1 / Rc;

    // devide into three cases
    if (D >= (PI * Rc / Bfe)) {
      sigma = sqrt(Bn / (2.0 * CN0) * D * (1.0 + 2.0 / (T * CN0 * (2 - D))));
    } else if (D > (Rc / Bfe)) {
      double tmp1 = Bn / (2.0 * CN0);
      double tmp2 = 1.0 / (Bfe * Tc) + Bfe * Tc / (PI - 1) * pow((D - 1.0 / (Bfe * Tc)), 2);
      double tmp3 = 1.0 + 2.0 / (T * CN0 * (2 - D));
      sigma = sqrt(tmp1 * tmp2 * tmp3);
    } else {
      sigma = sqrt(Bn / (2.0 * CN0) * (1.0 / (Bfe * Tc)) * (1.0 + 1.0 / (T * CN0)));
    }

    sigma = sigma * (C * Tc);  // convert to meters

    return sigma;
  }

  double GnssMeasurement::ComputeGnssPseudorangerateNoise(double CN0_dB, double lambda) {
    double F = 2;  // F=1 at high CN0, F=2 at low CN0
    double CN0 = pow(10, CN0_dB / 10);

    // extract gnss receiver parameters
    double Bn = gnssr_param.Bn;
    double T = gnssr_param.T;

    // compute pseudorange rate noise
    double sigma = lambda / (2 * PI * T) * sqrt(4 * F * Bn / CN0 * (1 + 1.0 / (T * CN0)));

    return sigma;
  }

  double GnssMeasurement::ComputeGnssCarrierPhaseNoise(double CN0_dB, double lambda) {
    // thermal noise in PLL
    double CN0 = pow(10, CN0_dB / 10);

    // extract gnss receiver parameters
    double Bp = gnssr_param.Bp;
    double T = gnssr_param.T;

    double sigma = lambda / (2 * PI) * sqrt(Bp / CN0 * (1.0 + 1.0 / (2 * T * CN0)));

    return sigma;
  }

}  // namespace lupnt
