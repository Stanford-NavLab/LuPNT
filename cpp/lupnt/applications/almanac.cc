#include "lupnt/applications/almanac.h"

#include <cmath>

#include "lupnt/applications/ephemeris_basis.h"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/core/error.h"

namespace lupnt {

  namespace {
    constexpr int kNumBaseParams = 2;  // t_ref, t_fit
    constexpr int kNumElements = 6;    // a, e, i, raan, argp, M

    int SegmentLen(int order) { return order + 1; }

    // Osculating classical orbital elements [a, e, i, raan, argp, M] at every
    // sample of `rv` (each row converted independently via `CartToClassical`),
    // with the angular elements unwrapped so they are safe to fit with a
    // polynomial.
    MatXd OsculatingElements(const MatXd& rv, double gm) {
      const int n = static_cast<int>(rv.rows());
      MatXd coe(n, 6);
      for (int i = 0; i < n; ++i) {
        const VecXd rv_row_d = rv.row(i).transpose();
        const Vec6 rv_row = rv_row_d.cast<Real>();
        const Vec6 coe_i = CartToClassical(rv_row, Real(gm));
        coe.row(i) = coe_i.cast<double>().transpose();
      }
      coe.col(2) = UnwrapAngles(coe.col(2));
      coe.col(3) = UnwrapAngles(coe.col(3));
      coe.col(4) = UnwrapAngles(coe.col(4));
      coe.col(5) = UnwrapAngles(coe.col(5));
      return coe;
    }
  }  // namespace

  VecXd Almanac::Fit(const VecXd& t_s, const MatXd& rv) const {
    const int n = static_cast<int>(t_s.size());
    const int order = options_.poly_order;
    LUPNT_CHECK(n >= order + 2, "Not enough samples to fit the requested polynomial order",
                "Almanac");
    LUPNT_CHECK(rv.rows() == n && rv.cols() == 6, "rv must be [N x 6]", "Almanac");

    const double t_ref = t_s(0);
    const double t_fit = t_s(n - 1) - t_s(0);
    LUPNT_CHECK(t_fit > 0.0, "t_s must be strictly increasing", "Almanac");
    const VecXd t_k = t_s.array() - t_ref;

    const MatXd coe = OsculatingElements(rv, options_.gm);
    const MatXd T = ChebyshevBasis(t_k, t_fit, order);
    auto lstsq = [&](const VecXd& y) -> VecXd { return T.colPivHouseholderQr().solve(y); };

    const VecXd a_coeffs = lstsq(coe.col(0));
    const VecXd e_coeffs = lstsq(coe.col(1));
    const VecXd i_coeffs = lstsq(coe.col(2));
    const VecXd raan_coeffs = lstsq(coe.col(3));
    const VecXd argp_coeffs = lstsq(coe.col(4));

    // Mean anomaly: remove the nominal two-body drift (integrated from the fitted
    // semi-major-axis history) before fitting a low-order polynomial residual.
    const VecXd a_series = T * a_coeffs;
    const VecXd n_series = (options_.gm * a_series.array().pow(-3)).sqrt();
    const VecXd m_nom = CumTrapz(n_series, t_k);
    const VecXd m_coeffs = lstsq(coe.col(5) - m_nom);

    const int seg = SegmentLen(order);
    VecXd params(kNumBaseParams + kNumElements * seg);
    params(0) = t_ref;
    params(1) = t_fit;
    params.segment(kNumBaseParams + 0 * seg, seg) = a_coeffs;
    params.segment(kNumBaseParams + 1 * seg, seg) = e_coeffs;
    params.segment(kNumBaseParams + 2 * seg, seg) = i_coeffs;
    params.segment(kNumBaseParams + 3 * seg, seg) = raan_coeffs;
    params.segment(kNumBaseParams + 4 * seg, seg) = argp_coeffs;
    params.segment(kNumBaseParams + 5 * seg, seg) = m_coeffs;
    return params;
  }

  MatXd Almanac::Eval(const VecXd& t_s, const VecXd& params) const {
    const int n = static_cast<int>(t_s.size());
    const int order = options_.poly_order;
    const int seg = SegmentLen(order);
    LUPNT_CHECK(params.size() == kNumBaseParams + kNumElements * seg,
                "params has the wrong size for this Almanac's options", "Almanac");

    const double t_ref = params(0);
    const double t_fit = params(1);
    const VecXd t_k = t_s.array() - t_ref;

    const MatXd T = ChebyshevBasis(t_k, t_fit, order);
    const VecXd a = T * params.segment(kNumBaseParams + 0 * seg, seg);
    const VecXd e = T * params.segment(kNumBaseParams + 1 * seg, seg);
    const VecXd inc = T * params.segment(kNumBaseParams + 2 * seg, seg);
    const VecXd raan = T * params.segment(kNumBaseParams + 3 * seg, seg);
    const VecXd argp = T * params.segment(kNumBaseParams + 4 * seg, seg);
    const VecXd m_res = T * params.segment(kNumBaseParams + 5 * seg, seg);

    const VecXd n_series = (options_.gm * a.array().pow(-3)).sqrt();
    const VecXd m_nom = CumTrapz(n_series, t_k);
    const VecXd m = m_nom + m_res;

    MatXd out(n, 6);
    for (int i = 0; i < n; ++i) {
      Vec6 coe_i;
      coe_i << a(i), e(i), inc(i), raan(i), argp(i), m(i);
      const Vec6 rv_i = ClassicalToCart(coe_i, Real(options_.gm));
      out.row(i) = rv_i.cast<double>().transpose();
    }
    return out;
  }

  EphemerisFitErrorStats Almanac::EvalError(const VecXd& t_s, const MatXd& rv_ref,
                                             const VecXd& params) const {
    return ComputeFitErrorStats(Eval(t_s, params), rv_ref);
  }

  int Almanac::NumParams() const { return kNumBaseParams + kNumElements * SegmentLen(options_.poly_order); }

  std::vector<std::string> Almanac::ParamNames() const {
    std::vector<std::string> names = {"t_ref", "t_fit"};
    for (const std::string elem : {"a", "e", "i", "raan", "argp", "M"}) {
      for (int j = 0; j <= options_.poly_order; ++j) names.push_back(elem + "_" + std::to_string(j));
    }
    return names;
  }

}  // namespace lupnt
