#include "lupnt/applications/ephemeris/lunanet_almanac.h"

#include <cmath>

#include "lupnt/applications/ephemeris/ephemeris_basis.h"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/core/error.h"

namespace lupnt {

  // LansAlmanac representation following Algorithm 2 of Iiyama & Gao, "Ephemeris and
  // LansAlmanac Design for Lunar Navigation Satellites": each osculating element is a
  // low-order polynomial plus an element-specific Fourier term, the argument of
  // periapsis is replaced by an argument-of-latitude correction u = nu + (poly +
  // Fourier) so that the reconstruction stays well-conditioned near periapsis of
  // eccentric orbits, and position/velocity are formed directly from
  // (a, e, i, node, u) with velocity by finite differencing.
  namespace {
    constexpr int kNumBaseParams = 3;  // t_ref, t_fit, a_ref
    constexpr int kNumElements = 6;    // a, e, i, node, M, u

    int SegmentLen(int poly_order, int num_fourier) { return (poly_order + 1) + 2 * num_fourier; }

    double FrameSpinRate(Frame frame) { return frame == Frame::MOON_PA ? OMEGA_MOON : 0.0; }
    Vec3d SpinCross(double spin, const Vec3d& r) { return spin * Vec3d(-r(1), r(0), 0.0); }

    // Base Fourier angular frequency [rad/s] for element `d`: the orbital mean
    // motion (2*pi/T_orb) for the semi-major axis (d==0), else 4*pi/T_sid.
    double ElementFreq(int d, double a_ref, double gm, double sidereal_period_s) {
      if (d == 0) return std::sqrt(gm / (a_ref * a_ref * a_ref));
      return 4.0 * PI / sidereal_period_s;
    }

    // Design matrix [N x SegmentLen]: 1, s, ..., s^p (s = t_k/t_fit) followed by
    // cos(h*freq*t_k), sin(h*freq*t_k) for h = 1..num_fourier.
    MatXd ElementBasis(const VecXd& t_k, double t_fit, int poly_order, int num_fourier,
                       double freq) {
      const int n = static_cast<int>(t_k.size());
      MatXd A(n, SegmentLen(poly_order, num_fourier));
      const VecXd s = t_k.array() / t_fit;
      A.col(0).setOnes();
      for (int p = 1; p <= poly_order; ++p) A.col(p) = A.col(p - 1).array() * s.array();
      const int off = poly_order + 1;
      for (int h = 1; h <= num_fourier; ++h) {
        A.col(off + 2 * (h - 1)) = (h * freq * t_k.array()).cos();
        A.col(off + 2 * (h - 1) + 1) = (h * freq * t_k.array()).sin();
      }
      return A;
    }

    double SolveKepler(double m, double e) {
      const double m_wrapped = std::atan2(std::sin(m), std::cos(m));
      double ecc_anom = m_wrapped;
      for (int it = 0; it < 60; ++it) {
        const double step
            = (ecc_anom - e * std::sin(ecc_anom) - m_wrapped) / (1.0 - e * std::cos(ecc_anom));
        ecc_anom -= step;
        if (std::abs(step) < 1e-13) break;
      }
      return ecc_anom;
    }

    double TrueAnomaly(double ecc_anom, double e) {
      return 2.0
             * std::atan2(std::sqrt(1.0 + e) * std::sin(0.5 * ecc_anom),
                          std::sqrt(1.0 - e) * std::cos(0.5 * ecc_anom));
    }

    // Osculating classical elements [a, e, i, node, argp, M] per sample, angles
    // unwrapped. If `spin` != 0 the states are in a frame rotating at `spin` about
    // z; the velocity is mapped to the inertial (PAI) frame (+omega x r) first.
    MatXd OsculatingElements(const MatXd& rv, double gm, double spin) {
      const int n = static_cast<int>(rv.rows());
      MatXd coe(n, 6);
      for (int i = 0; i < n; ++i) {
        Vec6d rv_row = rv.row(i).transpose();
        rv_row.tail<3>() += SpinCross(spin, rv_row.head<3>());
        const Vec6 rv_row_r = rv_row.cast<Real>();
        const Vec6 coe_i = CartToClassical(rv_row_r, Real(gm));
        coe.row(i) = coe_i.cast<double>().transpose();
      }
      coe.col(2) = UnwrapAngles(coe.col(2));  // i
      coe.col(3) = UnwrapAngles(coe.col(3));  // node
      coe.col(4) = UnwrapAngles(coe.col(4));  // argp
      coe.col(5) = UnwrapAngles(coe.col(5));  // M
      return coe;
    }

    // The six reconstructed element series (a, e, i, node, mean anomaly, and the
    // argument-of-latitude correction u_corr = argp) at `t_s` from the parameters.
    struct AlmanacSeries {
      VecXd a, e, inc, node, m, u_corr;
    };

    AlmanacSeries EvalSeries(const VecXd& t_s, const VecXd& params, const AlmanacFitOptions& opt) {
      const int poly = opt.poly_order;
      const int nf = std::max(0, opt.num_fourier_terms);
      const int seg = SegmentLen(poly, nf);
      const double t_ref = params(0), t_fit = params(1), a_ref = params(2);
      const VecXd t_k = t_s.array() - t_ref;
      const MatXd A_a = ElementBasis(t_k, t_fit, poly, nf,
                                     ElementFreq(0, a_ref, opt.gm, opt.sidereal_period_s));
      const MatXd A_s = ElementBasis(t_k, t_fit, poly, nf,
                                     ElementFreq(1, a_ref, opt.gm, opt.sidereal_period_s));
      AlmanacSeries s;
      s.a = A_a * params.segment(kNumBaseParams + 0 * seg, seg);
      s.e = A_s * params.segment(kNumBaseParams + 1 * seg, seg);
      s.inc = A_s * params.segment(kNumBaseParams + 2 * seg, seg);
      s.node = A_s * params.segment(kNumBaseParams + 3 * seg, seg);
      const VecXd m_res = A_s * params.segment(kNumBaseParams + 4 * seg, seg);
      s.u_corr = A_s * params.segment(kNumBaseParams + 5 * seg, seg);
      const VecXd n_series = (opt.gm * s.a.array().pow(-3)).sqrt();
      s.m = CumTrapz(n_series, t_k) + m_res;
      return s;
    }
  }  // namespace

  VecXd LansAlmanac::Fit(const VecXd& t_s, const MatXd& rv) const {
    const int n = static_cast<int>(t_s.size());
    const int poly = options_.poly_order;
    const int nf = std::max(0, options_.num_fourier_terms);
    const int seg = SegmentLen(poly, nf);
    LUPNT_CHECK(n >= seg, "Not enough samples to fit the requested polynomial/Fourier basis",
                "LansAlmanac");
    LUPNT_CHECK(rv.rows() == n && rv.cols() == 6, "rv must be [N x 6]", "LansAlmanac");

    const double t_ref = t_s(0);
    const double t_fit = t_s(n - 1) - t_s(0);
    LUPNT_CHECK(t_fit > 0.0, "t_s must be strictly increasing", "LansAlmanac");
    const VecXd t_k = t_s.array() - t_ref;

    const double spin = FrameSpinRate(options_.frame);
    const MatXd coe = OsculatingElements(rv, options_.gm, spin);
    const double a_ref = coe.col(0).mean();

    const MatXd A_a = ElementBasis(t_k, t_fit, poly, nf,
                                   ElementFreq(0, a_ref, options_.gm, options_.sidereal_period_s));
    const MatXd A_s = ElementBasis(t_k, t_fit, poly, nf,
                                   ElementFreq(1, a_ref, options_.gm, options_.sidereal_period_s));
    const auto qr_a = A_a.colPivHouseholderQr();
    const auto qr_s = A_s.colPivHouseholderQr();

    const VecXd a_coeffs = qr_a.solve(coe.col(0));
    const VecXd e_coeffs = qr_s.solve(coe.col(1));
    const VecXd i_coeffs = qr_s.solve(coe.col(2));
    const VecXd node_coeffs = qr_s.solve(coe.col(3));  // node fit directly (drifts in PA)

    const VecXd a_series = A_a * a_coeffs;
    const VecXd e_series = A_s * e_coeffs;
    const VecXd n_series = (options_.gm * a_series.array().pow(-3)).sqrt();
    const VecXd m_nom = CumTrapz(n_series, t_k);
    const VecXd m_coeffs = qr_s.solve(coe.col(5) - m_nom);
    const VecXd m_fit = m_nom + A_s * m_coeffs;

    // Argument-of-latitude correction u - nu(from fitted M), formed smoothly as
    // argp + wrap(nu_osc - nu_fit) so the least-squares target has no 2*pi jumps.
    VecXd u_target(n);
    for (int i = 0; i < n; ++i) {
      const double nu_o = TrueAnomaly(SolveKepler(coe(i, 5), coe(i, 1)), coe(i, 1));
      const double nu_f = TrueAnomaly(SolveKepler(m_fit(i), e_series(i)), e_series(i));
      const double dnu = std::atan2(std::sin(nu_o - nu_f), std::cos(nu_o - nu_f));
      u_target(i) = coe(i, 4) + dnu;  // argp (unwrapped) + small residual
    }
    const VecXd u_coeffs = qr_s.solve(u_target);

    VecXd params(kNumBaseParams + kNumElements * seg);
    params(0) = t_ref;
    params(1) = t_fit;
    params(2) = a_ref;
    params.segment(kNumBaseParams + 0 * seg, seg) = a_coeffs;
    params.segment(kNumBaseParams + 1 * seg, seg) = e_coeffs;
    params.segment(kNumBaseParams + 2 * seg, seg) = i_coeffs;
    params.segment(kNumBaseParams + 3 * seg, seg) = node_coeffs;
    params.segment(kNumBaseParams + 4 * seg, seg) = m_coeffs;
    params.segment(kNumBaseParams + 5 * seg, seg) = u_coeffs;
    return params;
  }

  MatXd LansAlmanac::Eval(const VecXd& t_s, const VecXd& params) const {
    const int n = static_cast<int>(t_s.size());
    const int seg = SegmentLen(options_.poly_order, std::max(0, options_.num_fourier_terms));
    LUPNT_CHECK(params.size() == kNumBaseParams + kNumElements * seg,
                "params has the wrong size for this LansAlmanac's options", "LansAlmanac");

    // Reconstruct the Cartesian state from (a, e, i, node, argp = u_corr, M): the
    // argument-of-latitude correction u_corr is exactly the effective argument of
    // periapsis (u = nu + u_corr), so ClassicalToCart with argp = u_corr yields
    // the same position as Algorithm 2's u-based reconstruction, plus the analytic
    // two-body velocity. For a rotating output frame, remove the omega x r offset.
    const AlmanacSeries s = EvalSeries(t_s, params, options_);
    const double spin = FrameSpinRate(options_.frame);
    MatXd out(n, 6);
    for (int i = 0; i < n; ++i) {
      // Wrap the (large, unwrapped) mean anomaly to [-pi, pi] so ClassicalToCart's
      // Kepler solver stays well-conditioned near periapsis of eccentric orbits.
      const double m_wrapped = std::atan2(std::sin(s.m(i)), std::cos(s.m(i)));
      Vec6 coe_i;
      coe_i << s.a(i), s.e(i), s.inc(i), s.node(i), s.u_corr(i), m_wrapped;
      const Vec6 rv_i = ClassicalToCart(coe_i, Real(options_.gm));
      Vec6d rv_out = rv_i.cast<double>();
      rv_out.tail<3>() -= SpinCross(spin, rv_out.head<3>());
      out.row(i) = rv_out.transpose();
    }
    return out;
  }

  EphemerisFitErrorStats LansAlmanac::EvalError(const VecXd& t_s, const MatXd& rv_ref,
                                                const VecXd& params) const {
    return ComputeFitErrorStats(Eval(t_s, params), rv_ref);
  }

  int LansAlmanac::NumParams() const {
    return kNumBaseParams
           + kNumElements
                 * SegmentLen(options_.poly_order, std::max(0, options_.num_fourier_terms));
  }

  std::vector<std::string> LansAlmanac::ParamNames() const {
    const int num_fourier = std::max(0, options_.num_fourier_terms);
    std::vector<std::string> names = {"t_ref", "t_fit", "a_ref"};
    for (const std::string elem : {"a", "e", "i", "raan", "M", "u"}) {
      for (int p = 0; p <= options_.poly_order; ++p)
        names.push_back(elem + "_p" + std::to_string(p));
      for (int h = 1; h <= num_fourier; ++h) {
        names.push_back(elem + "_fc" + std::to_string(h));
        names.push_back(elem + "_fs" + std::to_string(h));
      }
    }
    return names;
  }

}  // namespace lupnt
