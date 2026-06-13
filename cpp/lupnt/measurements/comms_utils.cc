#include "lupnt/measurements/comms_utils.h"

#include "lupnt/core/constants.h"

namespace lupnt {

  Real SigmaDll(const DllParams& params, Real CN0_w) {
    Real B_dll = params.B_dll;
    Real T_i = params.T_i;
    Real d = params.d;
    Real Tc = params.Tc;
    Real B_fe = params.B_fe;

    Real term_1 = B_dll / (2.0 * CN0_w);
    Real term_2 = 0.0;
    Real term_3 = 1.0 + 2.0 / (T_i * CN0_w * (2.0 - d));

    Real sigma_dll_2 = 0.0;
    if (d <= 1.0 / (Tc * B_fe)) {
      term_2 = 1.0 / (Tc * B_fe);
      term_3 = 1.0 + (1.0 / (T_i * CN0_w));
      sigma_dll_2 = term_1 * term_2 * term_3;
    } else if ((PI / (Tc * B_fe) > d) && (d > 1.0 / (Tc * B_fe))) {
      term_2 = 1.0 / (Tc * B_fe) + Tc * B_fe / (PI - 1) * pow(d - 1.0 / (Tc * B_fe), 2);
      sigma_dll_2 = term_1 * term_2 * term_3;
    } else {
      sigma_dll_2 = term_1 * d * term_3;
    }
    return sqrt(sigma_dll_2);
  }

  ArrX SigmaDll(const DllParams& params, const ArrX& CN0_w) {
    Real B_dll = params.B_dll;
    Real T_i = params.T_i;
    Real d = params.d;
    Real Tc = params.Tc;
    Real B_fe = params.B_fe;

    int n = CN0_w.rows();
    int m = CN0_w.cols();

    ArrX ones = ArrX::Ones(n, m);
    ArrX zeros = ArrX::Zero(n, m);

    ArrX term_1 = B_dll / (2.0 * CN0_w);
    ArrX term_2 = zeros;
    ArrX term_3 = 1.0 + 2.0 / (T_i * CN0_w * (2.0 - d));

    ArrX sigma_dll_2 = zeros;
    if (d <= 1.0 / (Tc * B_fe)) {
      term_2 = ones / (Tc * B_fe);
      term_3 = 1.0 + (1.0 / (T_i * CN0_w));
      sigma_dll_2 = term_1 * term_2 * term_3;
    } else if ((PI / (Tc * B_fe) > d) && (d > 1.0 / (Tc * B_fe))) {
      term_2 = ones / (Tc * B_fe) + Tc * B_fe / (PI - 1) * pow(d - 1.0 / (Tc * B_fe), 2);
      sigma_dll_2 = term_1 * term_2 * term_3;
    } else {
      sigma_dll_2 = term_1 * d * term_3;
    }
    return sigma_dll_2.sqrt();
  }

  Real SigmaPll(const PllParams& params, Real CN0_w) {
    Real B_pll = params.B_pll;
    Real T_i = params.T_i;

    Real sigma_pll_2 = B_pll / CN0_w * (1.0 + 1.0 / (2.0 * T_i * CN0_w));
    return sqrt(sigma_pll_2);
  }

  ArrX SigmaPll(const PllParams& params, const ArrX& CN0_w) {
    Real B_pll = params.B_pll;
    Real T_i = params.T_i;

    ArrX sigma_pll_2 = B_pll / CN0_w * (1.0 + 1.0 / (2.0 * T_i * CN0_w));
    return sigma_pll_2.sqrt();
  }

  Real SigmaFll(const FllParams& params, Real CN0_w) {
    Real B_fll = params.B_fll;
    Real T_i = params.T_i;

    Real CN0_F_fll_w = DecibelToDecimal(params.CN0_F_fll);
    Real F_ll = (CN0_w >= CN0_F_fll_w) ? 1.0 : 2.0;

    Real sigma_fll_2 = 4.0 * F_ll * B_fll * (1.0 + 1.0 / (T_i * CN0_w));
    return sqrt(sigma_fll_2);
  }

  ArrX SigmaFll(const FllParams& params, const ArrX& CN0_w) {
    Real B_fll = params.B_fll;
    Real T_i = params.T_i;

    ArrX ones = ArrX::Ones(CN0_w.rows(), CN0_w.cols());
    Real CN0_F_fll_w = DecibelToDecimal(params.CN0_F_fll);
    ArrX F_ll = (CN0_w >= CN0_F_fll_w).select(ones, 2.0 * ones);

    ArrX sigma_fll_2 = 4.0 * F_ll * B_fll * (1.0 + 1.0 / (T_i * CN0_w));
    return sigma_fll_2.sqrt();
  }

  Real FreeSpacePathLoss(Real dist, Real freq) {
    return 20.0 * log10((4.0 * PI * dist) / (C / freq));
  }

  ArrX FreeSpacePathLoss(const ArrX& dist, Real freq) {
    return 20.0 * ((4.0 * PI * dist) / (C / freq)).log10();
  }

  Real ParabolicAntennaGain(Real phi, Real hpbw, Real G_max) {
    return G_max - 12.0 * phi * phi / pow(hpbw, 2);
  }

  ArrX ParabolicAntennaGain(const ArrX& phi, Real hpbw, Real G_max) {
    return G_max - 12.0 * phi.square() / pow(hpbw, 2);
  }

}  // namespace lupnt
