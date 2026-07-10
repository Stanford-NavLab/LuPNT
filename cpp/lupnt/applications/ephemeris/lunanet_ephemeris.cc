#include "lupnt/applications/ephemeris/lunanet_ephemeris.h"

#include <cmath>

#include "lupnt/applications/ephemeris/ephemeris_basis.h"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/core/error.h"

namespace lupnt {

  namespace {
    constexpr int kNumBaseParamsKepler = 8;    // t_ref, t_fit, a, e, i, raan, argp, M_ref
    constexpr int kNumBaseParamsNoKepler = 2;  // t_ref, t_fit

    int BaseParamCount(bool use_keplerian_baseline) {
      return use_keplerian_baseline ? kNumBaseParamsKepler : kNumBaseParamsNoKepler;
    }

    // Angular velocity [rad/s] of `frame` about its z-axis relative to inertial:
    // OMEGA_MOON for the rotating Moon-fixed principal-axis frame, 0 otherwise
    // (so an inertial `frame`, e.g. MOON_CI, reproduces a standard Kepler orbit).
    double FrameSpinRate(Frame frame) { return frame == Frame::MOON_PA ? OMEGA_MOON : 0.0; }

    // z-axis (Moon pole) cross product omega x r, with omega = spin * z_hat.
    Vec3d SpinCross(double spin, const Vec3d& r) { return spin * Vec3d(-r(1), r(0), 0.0); }

    // Result of propagating the osculating two-body baseline: Cartesian state in
    // the (possibly rotating) output frame, plus the argument of latitude u and
    // its rate (used only by the optional Fourier terms).
    struct BaselineArc {
      MatXd rv;    // [N x 6]
      VecXd u;     // [N] argument of latitude [rad]
      VecXd udot;  // [N] du/dt [rad/s]
    };

    double SolveKepler(double m, double e) {
      const double m_wrapped = std::fmod(m, 2.0 * PI);
      double ecc_anom = m_wrapped;
      for (int it = 0; it < 50; ++it) {
        const double f = ecc_anom - e * std::sin(ecc_anom) - m_wrapped;
        const double fp = 1.0 - e * std::cos(ecc_anom);
        const double step = f / fp;
        ecc_anom -= step;
        if (std::abs(step) < 1e-13) break;
      }
      return ecc_anom;
    }

    // Osculating two-body Kepler orbit [a, e, i, raan, argp, M_ref] fixed at
    // t_k=0, evaluated at every `t_k`. If `spin` != 0 the baseline is expressed in
    // a frame rotating about the z-axis at rate `spin`: the inertial two-body
    // motion is rigidly rotated by R_z(-spin*t_k), so the ascending node drifts at
    // -spin and the velocity carries the -spin x r term (Algorithm 1 of Iiyama &
    // Gao). Returns zeros (and empty u/udot) if `coe` is empty.
    BaselineArc EvalKeplerianBaseline(const VecXd& t_k, const VecXd& coe, double gm, double spin,
                                      bool want_arg_lat) {
      const int n = static_cast<int>(t_k.size());
      BaselineArc arc;
      arc.rv = MatXd::Zero(n, 6);
      if (coe.size() == 0) return arc;
      const double a = coe(0), e = coe(1), inc = coe(2), raan = coe(3), argp = coe(4),
                   m_ref = coe(5);
      const double mean_motion = std::sqrt(gm / (a * a * a));
      if (want_arg_lat) {
        arc.u = VecXd::Zero(n);
        arc.udot = VecXd::Zero(n);
      }
      for (int i = 0; i < n; ++i) {
        Vec6 coe_k;
        coe_k << a, e, inc, raan, argp, m_ref + mean_motion * t_k(i);
        const Vec6 rv_k = ClassicalToCart(coe_k, Real(gm));
        const Vec6d rv_pai = rv_k.cast<double>();

        // Rotate the inertial (PAI) two-body state into the rotating output frame.
        Vec6d rv_out = rv_pai;
        if (spin != 0.0) {
          const double th = -spin * t_k(i);
          const double c = std::cos(th), s = std::sin(th);
          Mat3d C;
          C << c, -s, 0, s, c, 0, 0, 0, 1;
          Mat3d Cdot;  // dC/dt = (dC/dtheta) * (-spin)
          Cdot << s, c, 0, -c, s, 0, 0, 0, 0;
          Cdot *= spin;
          const Vec3d r_pai = rv_pai.head<3>();
          const Vec3d v_pai = rv_pai.tail<3>();
          rv_out.head<3>() = C * r_pai;
          rv_out.tail<3>() = C * v_pai + Cdot * r_pai;
        }
        arc.rv.row(i) = rv_out.transpose();

        if (want_arg_lat) {
          const double ecc_anom = SolveKepler(m_ref + mean_motion * t_k(i), e);
          const double nu = 2.0
                            * std::atan2(std::sqrt(1.0 + e) * std::sin(0.5 * ecc_anom),
                                         std::sqrt(1.0 - e) * std::cos(0.5 * ecc_anom));
          const double edot = mean_motion / (1.0 - e * std::cos(ecc_anom));
          arc.u(i) = argp + nu;
          arc.udot(i) = edot * std::sqrt(1.0 - e * e) / (1.0 - e * std::cos(ecc_anom));
        }
      }
      return arc;
    }

    // Fourier design matrix [N x 2M]: columns [cos(u), sin(u), cos(2u), sin(2u), ...].
    MatXd FourierBasis(const VecXd& u, int num_terms) {
      const int n = static_cast<int>(u.size());
      MatXd F(n, 2 * num_terms);
      for (int h = 1; h <= num_terms; ++h) {
        F.col(2 * (h - 1)) = (h * u.array()).cos();
        F.col(2 * (h - 1) + 1) = (h * u.array()).sin();
      }
      return F;
    }

    // Time derivative of FourierBasis given u(t) and udot(t).
    MatXd FourierBasisDt(const VecXd& u, const VecXd& udot, int num_terms) {
      const int n = static_cast<int>(u.size());
      MatXd Fd(n, 2 * num_terms);
      for (int h = 1; h <= num_terms; ++h) {
        Fd.col(2 * (h - 1)) = (-h) * (h * u.array()).sin() * udot.array();
        Fd.col(2 * (h - 1) + 1) = h * (h * u.array()).cos() * udot.array();
      }
      return Fd;
    }
  }  // namespace

  VecXd CartesianEphemeris::Fit(const VecXd& t_s, const MatXd& rv) const {
    const int n = static_cast<int>(t_s.size());
    const int order = options_.poly_order;
    const int num_fourier = std::max(0, options_.num_fourier_terms);
    LUPNT_CHECK(num_fourier == 0 || options_.use_keplerian_baseline,
                "Fourier terms require use_keplerian_baseline", "CartesianEphemeris");
    LUPNT_CHECK(n >= order + 2 * num_fourier + 2,
                "Not enough samples to fit the requested Chebyshev/Fourier basis",
                "CartesianEphemeris");
    LUPNT_CHECK(rv.rows() == n && rv.cols() == 6, "rv must be [N x 6]", "CartesianEphemeris");

    const double spin = FrameSpinRate(options_.frame);
    const int idx_mid = n / 2;
    const double t_ref = t_s(idx_mid);
    const double t_fit = t_s(n - 1) - t_s(0);
    LUPNT_CHECK(t_fit > 0.0, "t_s must be strictly increasing", "CartesianEphemeris");
    const VecXd t_k = t_s.array() - t_ref;

    VecXd coe_ref;  // empty unless a Keplerian baseline is used
    if (options_.use_keplerian_baseline) {
      // Osculating elements from the reference-epoch state, in the inertial (PAI)
      // frame: undo the frame rotation's -omega x r velocity offset (position is
      // identical in the rotating and inertial frames at the reference epoch).
      Vec6d rv_ref = rv.row(idx_mid).transpose();
      rv_ref.tail<3>() += SpinCross(spin, rv_ref.head<3>());
      const Vec6 rv_ref_r = rv_ref.cast<Real>();
      const Vec6 coe_ref_r = CartToClassical(rv_ref_r, Real(options_.gm));
      coe_ref = coe_ref_r.cast<double>();
    }
    const BaselineArc arc = EvalKeplerianBaseline(t_k, coe_ref, options_.gm, spin, num_fourier > 0);
    const MatXd residual_pos = rv.leftCols(3) - arc.rv.leftCols(3);

    const int cheb_len = order + 1;
    const int four_len = 2 * num_fourier;
    const MatXd T = ChebyshevBasis(t_k, t_fit, order);
    const MatXd F = num_fourier > 0 ? FourierBasis(arc.u, num_fourier) : MatXd(n, 0);

    // Joint least-squares solve over [Chebyshev | Fourier] for each axis.
    MatXd design(n, cheb_len + four_len);
    design.leftCols(cheb_len) = T;
    if (four_len > 0) design.rightCols(four_len) = F;
    const auto qr = design.colPivHouseholderQr();
    MatXd cheb_coeffs(cheb_len, 3);
    MatXd four_coeffs(four_len, 3);
    for (int d = 0; d < 3; ++d) {
      const VecXd sol = qr.solve(residual_pos.col(d));
      cheb_coeffs.col(d) = sol.head(cheb_len);
      if (four_len > 0) four_coeffs.col(d) = sol.tail(four_len);
    }

    const int base = BaseParamCount(options_.use_keplerian_baseline);
    VecXd params(base + 3 * cheb_len + 3 * four_len);
    params(0) = t_ref;
    params(1) = t_fit;
    if (options_.use_keplerian_baseline) params.segment(2, 6) = coe_ref;
    for (int d = 0; d < 3; ++d) params.segment(base + d * cheb_len, cheb_len) = cheb_coeffs.col(d);
    for (int d = 0; d < 3; ++d) {
      if (four_len > 0)
        params.segment(base + 3 * cheb_len + d * four_len, four_len) = four_coeffs.col(d);
    }
    return params;
  }

  MatXd CartesianEphemeris::Eval(const VecXd& t_s, const VecXd& params) const {
    const int n = static_cast<int>(t_s.size());
    const int order = options_.poly_order;
    const int num_fourier = std::max(0, options_.num_fourier_terms);
    const int base = BaseParamCount(options_.use_keplerian_baseline);
    const int cheb_len = order + 1;
    const int four_len = 2 * num_fourier;
    LUPNT_CHECK(params.size() == base + 3 * cheb_len + 3 * four_len,
                "params has the wrong size for this CartesianEphemeris's options",
                "CartesianEphemeris");

    const double spin = FrameSpinRate(options_.frame);
    const double t_ref = params(0);
    const double t_fit = params(1);
    const VecXd t_k = t_s.array() - t_ref;

    const VecXd coe_ref = options_.use_keplerian_baseline ? VecXd(params.segment(2, 6)) : VecXd();
    const BaselineArc arc = EvalKeplerianBaseline(t_k, coe_ref, options_.gm, spin, num_fourier > 0);

    const MatXd T = ChebyshevBasis(t_k, t_fit, order);
    const MatXd Tdot = ChebyshevBasisDt(t_k, t_fit, order);
    const MatXd F = num_fourier > 0 ? FourierBasis(arc.u, num_fourier) : MatXd(n, 0);
    const MatXd Fdot = num_fourier > 0 ? FourierBasisDt(arc.u, arc.udot, num_fourier) : MatXd(n, 0);

    MatXd out(n, 6);
    for (int d = 0; d < 3; ++d) {
      const VecXd c = params.segment(base + d * cheb_len, cheb_len);
      out.col(d) = arc.rv.col(d) + T * c;
      out.col(3 + d) = arc.rv.col(3 + d) + Tdot * c;
      if (four_len > 0) {
        const VecXd cf = params.segment(base + 3 * cheb_len + d * four_len, four_len);
        out.col(d) += F * cf;
        out.col(3 + d) += Fdot * cf;
      }
    }
    return out;
  }

  EphemerisFitErrorStats CartesianEphemeris::EvalError(const VecXd& t_s, const MatXd& rv_ref,
                                                       const VecXd& params) const {
    return ComputeFitErrorStats(Eval(t_s, params), rv_ref);
  }

  int CartesianEphemeris::NumParams() const {
    const int num_fourier = std::max(0, options_.num_fourier_terms);
    return BaseParamCount(options_.use_keplerian_baseline) + 3 * (options_.poly_order + 1)
           + 3 * 2 * num_fourier;
  }

  std::vector<std::string> CartesianEphemeris::ParamNames() const {
    std::vector<std::string> names = {"t_ref", "t_fit"};
    if (options_.use_keplerian_baseline) {
      for (const std::string s : {"a", "e", "i", "raan", "argp", "M_ref"}) names.push_back(s);
    }
    for (const std::string axis : {"x", "y", "z"}) {
      for (int j = 0; j <= options_.poly_order; ++j)
        names.push_back(axis + "_" + std::to_string(j));
    }
    const int num_fourier = std::max(0, options_.num_fourier_terms);
    for (const std::string axis : {"x", "y", "z"}) {
      for (int h = 1; h <= num_fourier; ++h) {
        names.push_back(axis + "_fc_" + std::to_string(h));
        names.push_back(axis + "_fs_" + std::to_string(h));
      }
    }
    return names;
  }

}  // namespace lupnt
