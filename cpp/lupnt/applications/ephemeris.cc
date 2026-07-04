#include "lupnt/applications/ephemeris.h"

#include <cmath>

#include "lupnt/applications/ephemeris_basis.h"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/core/error.h"

namespace lupnt {

  namespace {
    constexpr int kNumBaseParamsKepler = 8;    // t_ref, t_fit, a, e, i, raan, argp, M_ref
    constexpr int kNumBaseParamsNoKepler = 2;  // t_ref, t_fit

    int BaseParamCount(bool use_keplerian_baseline) {
      return use_keplerian_baseline ? kNumBaseParamsKepler : kNumBaseParamsNoKepler;
    }

    // Osculating two-body Kepler orbit [a, e, i, raan, argp, M_ref] fixed at
    // t_k=0, evaluated at every `t_k` via `ClassicalToCart`. Returns an [N x 6]
    // matrix of zeros if `coe` is empty (i.e. no Keplerian baseline requested).
    MatXd EvalKeplerianBaseline(const VecXd& t_k, const VecXd& coe, double gm) {
      const int n = static_cast<int>(t_k.size());
      MatXd rv_baseline = MatXd::Zero(n, 6);
      if (coe.size() == 0) return rv_baseline;
      const double a = coe(0), e = coe(1), inc = coe(2), raan = coe(3), argp = coe(4),
                   m_ref = coe(5);
      const double mean_motion = std::sqrt(gm / (a * a * a));
      for (int i = 0; i < n; ++i) {
        Vec6 coe_k;
        coe_k << a, e, inc, raan, argp, m_ref + mean_motion * t_k(i);
        const Vec6 rv_k = ClassicalToCart(coe_k, Real(gm));
        rv_baseline.row(i) = rv_k.cast<double>().transpose();
      }
      return rv_baseline;
    }
  }  // namespace

  VecXd CartesianEphemeris::Fit(const VecXd& t_s, const MatXd& rv) const {
    const int n = static_cast<int>(t_s.size());
    const int order = options_.poly_order;
    LUPNT_CHECK(n >= order + 2, "Not enough samples to fit the requested Chebyshev order",
                "CartesianEphemeris");
    LUPNT_CHECK(rv.rows() == n && rv.cols() == 6, "rv must be [N x 6]", "CartesianEphemeris");

    const int idx_mid = n / 2;
    const double t_ref = t_s(idx_mid);
    const double t_fit = t_s(n - 1) - t_s(0);
    LUPNT_CHECK(t_fit > 0.0, "t_s must be strictly increasing", "CartesianEphemeris");
    const VecXd t_k = t_s.array() - t_ref;

    VecXd coe_ref;  // empty unless a Keplerian baseline is used
    if (options_.use_keplerian_baseline) {
      const VecXd rv_ref_d = rv.row(idx_mid).transpose();
      const Vec6 rv_ref = rv_ref_d.cast<Real>();
      const Vec6 coe_ref_r = CartToClassical(rv_ref, Real(options_.gm));
      coe_ref = coe_ref_r.cast<double>();
    }
    const MatXd rv_baseline = EvalKeplerianBaseline(t_k, coe_ref, options_.gm);
    const MatXd residual_pos = rv.leftCols(3) - rv_baseline.leftCols(3);

    const MatXd T = ChebyshevBasis(t_k, t_fit, order);
    MatXd coeffs(order + 1, 3);
    for (int d = 0; d < 3; ++d) coeffs.col(d) = T.colPivHouseholderQr().solve(residual_pos.col(d));

    const int base = BaseParamCount(options_.use_keplerian_baseline);
    VecXd params(base + 3 * (order + 1));
    params(0) = t_ref;
    params(1) = t_fit;
    if (options_.use_keplerian_baseline) params.segment(2, 6) = coe_ref;
    for (int d = 0; d < 3; ++d) params.segment(base + d * (order + 1), order + 1) = coeffs.col(d);
    return params;
  }

  MatXd CartesianEphemeris::Eval(const VecXd& t_s, const VecXd& params) const {
    const int n = static_cast<int>(t_s.size());
    const int order = options_.poly_order;
    const int base = BaseParamCount(options_.use_keplerian_baseline);
    LUPNT_CHECK(params.size() == base + 3 * (order + 1),
                "params has the wrong size for this CartesianEphemeris's options",
                "CartesianEphemeris");

    const double t_ref = params(0);
    const double t_fit = params(1);
    const VecXd t_k = t_s.array() - t_ref;

    const VecXd coe_ref
        = options_.use_keplerian_baseline ? VecXd(params.segment(2, 6)) : VecXd();
    const MatXd rv_baseline = EvalKeplerianBaseline(t_k, coe_ref, options_.gm);

    const MatXd T = ChebyshevBasis(t_k, t_fit, order);
    const MatXd Tdot = ChebyshevBasisDt(t_k, t_fit, order);

    MatXd out(n, 6);
    for (int d = 0; d < 3; ++d) {
      const VecXd c = params.segment(base + d * (order + 1), order + 1);
      out.col(d) = rv_baseline.col(d) + T * c;
      out.col(3 + d) = rv_baseline.col(3 + d) + Tdot * c;
    }
    return out;
  }

  EphemerisFitErrorStats CartesianEphemeris::EvalError(const VecXd& t_s, const MatXd& rv_ref,
                                                        const VecXd& params) const {
    return ComputeFitErrorStats(Eval(t_s, params), rv_ref);
  }

  int CartesianEphemeris::NumParams() const {
    return BaseParamCount(options_.use_keplerian_baseline) + 3 * (options_.poly_order + 1);
  }

  std::vector<std::string> CartesianEphemeris::ParamNames() const {
    std::vector<std::string> names = {"t_ref", "t_fit"};
    if (options_.use_keplerian_baseline) {
      for (const std::string s : {"a", "e", "i", "raan", "argp", "M_ref"}) names.push_back(s);
    }
    for (const std::string axis : {"x", "y", "z"}) {
      for (int j = 0; j <= options_.poly_order; ++j) names.push_back(axis + "_" + std::to_string(j));
    }
    return names;
  }

}  // namespace lupnt
