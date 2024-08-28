
#include "lupnt/measurements/comm_utils.h"

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {
  double ComputeEsN0(double EbN0, Modulation modulation_type, double coding_rate) {
    double K = BitsPerSymbol(modulation_type);
    double EsN0 = EbN0 * K * coding_rate;
    return EsN0;
  }

  double BitsPerSymbol(Modulation modulation_type) {
    double M = 1.0;
    switch (modulation_type) {
      case Modulation::Residual:
      case Modulation::BPSK:
      case Modulation::GMSK:
      case Modulation::GMSK_PN: M = 1; break;
      case Modulation::QPSK:
      case Modulation::OQPSK: M = 2; break;
      default: std::cerr << "Modulation type not supported" << std::endl;
    }
    return M;
  }

  double ComputeBER(double EbN0, Modulation modulation_type) {
    double BER = 0.0;

    switch (modulation_type) {
      case Modulation::BPSK:
      case Modulation::QPSK:
      case Modulation::OQPSK:
      case Modulation::GMSK:
      case Modulation::GMSK_PN: BER = 0.5 * erfc(sqrt(EbN0)); break;

      default: std::cerr << "Modulation type not supported" << std::endl;
    }

    return BER;
  }

  FrequencyBand GetFrequencyBand(double f_C) {
    FrequencyBand fbu = FrequencyBand::S;

    // UHF, S, X, Ku, Ka
    if (f_C >= 300e6 && f_C < 1e9) {
      fbu = FrequencyBand::UHF;
    } else if (f_C >= 1.0e9 && f_C < 2.0e9) {
      fbu = FrequencyBand::L;
    } else if (f_C >= 2e9 && f_C < 4e9) {
      fbu = FrequencyBand::S;
    } else if (f_C >= 2e9 && f_C < 4e9) {
      fbu = FrequencyBand::Cband;
    } else if (f_C >= 8e9 && f_C < 12e9) {
      fbu = FrequencyBand::X;
    } else if (f_C >= 12e9 && f_C < 18e9) {
      fbu = FrequencyBand::Ku;
    } else if (f_C >= 18e9 && f_C < 26e9) {
      fbu = FrequencyBand::K;
    } else if (f_C >= 26e9 && f_C < 40e9) {
      fbu = FrequencyBand::Ka;
    }

    return fbu;
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
                               Modulation modulation_type, double m_R) {
    double S_L = 1.0;
    double EsN0 = PT_N0 * T_s;
    double rho_L = PT_N0 / B_L_carrier;
    double tmp1, tmp2, tmp3;

    // Compute carier loop signal-to-noise ratio
    if (modulation_type == Modulation::Residual) {
      // carrier nominal power https://public.ccsds.org/Pubs/401x0b17s.pdf
      double PC_N0 = PT_N0 * cos(m_R) * pow(J0Bessel(m_R), 2);
      rho_L = PC_N0 / B_L_carrier * EsN0 / (1 + 2 * EsN0);

    } else {
      switch (modulation_type) {
        case Modulation::BPSK: S_L = 2 * EsN0 / (1 + 2 * EsN0); break;
        case Modulation::QPSK:
          tmp1 = (9 / 4) / EsN0;
          tmp2 = (3 / 2) / pow(EsN0, 2);
          tmp3 = (3 / 16) / pow(EsN0, 3);

          S_L = 1 / (1 + tmp1 + tmp2 + tmp3);
          break;
        case Modulation::OQPSK:
        case Modulation::GMSK:
        case Modulation::GMSK_PN:
          // same as OQPSK
          // Reference: GMSK Modulation for Deep Space Application
          // https://ieeexplore-ieee-org.stanford.idm.oclc.org/stamp/stamp.jsp?tp=&arnumber=6187097
          tmp1 = (9 / 4) / EsN0;
          tmp2 = (3 / 2) / pow(EsN0, 2);
          tmp3 = (3 / 16) / pow(EsN0, 3);

          S_L = 1 / 4 / (1 + tmp1 + tmp2 + tmp3);
          break;
        case Modulation::Residual: throw std::runtime_error("Not implemented");
      }
      rho_L = PT_N0 * S_L / B_L_carrier;
    }

    // Recommended value to avoid cycle slips
    // rho_L >= 10 dB for residual carrier
    // rho_L >= 17 dB for BPSK
    // rho_L >= 23 dB for QPSK and Offset QPSK

    return rho_L;
  }

}  // namespace lupnt
