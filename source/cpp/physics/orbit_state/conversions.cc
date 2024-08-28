#include "lupnt/physics/orbit_state/conversions.h"

#include <cmath>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state/anomaly.h"
#include "lupnt/physics/orbit_state/orbit_states.h"

namespace lupnt {

  /// @brief Convert Cartesian state to classical orbital elements
  /// @param rv Cartesian state [km, km/s]
  /// @param GM Gravitational parameter [km^3/s^2]
  /// @return Classical orbital elements [km, –, rad, rad, rad, rad]
  ClassicalOE Cart2Classical(const CartesianOrbitState &rv, Real GM) {
    return ClassicalOE(Cart2Classical(rv.GetVec6(), GM), rv.GetFrame());
  }

  /// @brief Convert Cartesian state to classical orbital elements
  /// @param rv Cartesian state [km, km/s]
  /// @param GM Gravitational parameter [km^3/s^2]
  /// @return Classical orbital elements [km, –, rad, rad, rad, rad]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec6 Cart2Classical(const Vec6 &rv, Real GM) {
    Vec3 r = rv.head(3);  // Position
    Vec3 v = rv.tail(3);  // Velocity
    Vec3 h = r.cross(v);  // Areal velocity
    Real H = h.norm();    // Angular momentum

    Real Omega = atan2(h(0), -h(1));                        // Long. ascend. node
    Real i = atan2(sqrt(h(0) * h(0) + h(1) * h(1)), h(2));  // Inclination
    Real u = atan2(r(2) * H, -r(0) * h(1) + r(1) * h(0));   // Arg. of latitude
    Real R = r.norm();                                      // Distance
    Real a = 1.0 / (2.0 / R - v.squaredNorm() / GM);        // Semi-major axis

    Real eCosE = 1.0 - R / a;              // e*cos(E)
    Real eSinE = r.dot(v) / sqrt(GM * a);  // e*sin(E)
    Real e2 = eCosE * eCosE + eSinE * eSinE;
    Real e = sqrt(e2);             // Eccentricity
    Real E = atan2(eSinE, eCosE);  // Eccentric anomaly

    Real M = Wrap2Pi(E - eSinE);                          // Mean anomaly
    Real nu = atan2(sqrt(1.0 - e2) * eSinE, eCosE - e2);  // True anomaly
    Real omega = Wrap2Pi(u - nu);                         // Arg. of perihelion

    return Vec6(a, e, i, Omega, omega, M);
  }

  /// @brief Compute the classical orbital elements from two position vectors and
  /// the time difference between them
  /// @param dt Time difference between the two position vectors [s]
  /// @param r1 Initial position vector [km]
  /// @param r2 Final position vector [km]
  /// @param GM Gravitational parameter [km^3/s^2]
  /// @return Classical orbital elements [km, –, rad, rad, rad, rad]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec6 Cart2Classical(Real dt, const Vec3 &r1, const Vec3 &r2, Real GM) {
    // Calculate vector r0 (fraction of r2 perpendicular to r1)
    // and the magnitudes of r1, r2 and r0
    Real s_a = r1.norm();
    Vec3 e1 = r1 / s_a;
    Real s_b = r2.norm();
    Real fac = r2.dot(e1);
    Vec3 r0 = r2 - fac * e1;
    Real s0 = r0.norm();
    Vec3 e0 = r0 / s0;

    // Inclination and ascending node
    Vec3 W = e1.cross(e0);
    Real Omega = (atan2(W(0), -W(1)));                      // Long. ascend. node
    Real i = atan2(sqrt(W(0) * W(0) + W(1) * W(1)), W(2));  // Inclination
    Real u;
    if (i == 0.0)
      u = atan2(r1(1), r1(0));
    else
      u = atan2(+e1(2), -e1(0) * W(1) + e1(1) * W(0));

    // Semilatus rectum
    Real tau = sqrt(GM) * dt;
    Real eta = RatioOfSectorToTriangleArea(r1, r2, tau);
    Real p = pow(s_a * s0 * eta / tau, 2);

    // Eccentricity, true anomaly and argument of perihelion
    Real cos_dnu = fac / s_b;
    Real sin_dnu = s0 / s_b;

    Real ecos_nu = p / s_a - 1.0;
    Real esin_nu = (ecos_nu * cos_dnu - (p / s_b - 1.0)) / sin_dnu;

    Real e = sqrt(ecos_nu * ecos_nu + esin_nu * esin_nu);
    Real nu = atan2(esin_nu, ecos_nu);
    Real omega = Wrap2Pi(u - nu);

    // Perihelion distance, semimajor axis and mean motion
    Real a = p / (1.0 - e * e);

    // Mean anomaly and time of perihelion passage
    Real M;
    if (e < 1.0) {
      Real E = atan2(sqrt((1.0 - e) * (1.0 + e)) * esin_nu, ecos_nu + e * e);
      M = Wrap2Pi(E - e * sin(E));
    } else {
      Real sinhH = sqrt((e - 1.0) * (e + 1.0)) * esin_nu / (e + e * ecos_nu);
      M = e * sinhH - log(sinhH + sqrt(1.0 + sinhH * sinhH));
    }

    // Keplerian elements vector
    return Vec6(a, e, i, Omega, omega, M);
  }

  /// @brief Convert two cartesian states to a relative RTN state
  /// @param rv_c Chief spacecraft state [km, km/s]
  /// @param rv_d Deputy spacecraft state [km, km/s]
  /// @return Relative RTN state [km, km/s]
  CartesianOrbitState Inertial2Synodic(const CartesianOrbitState &rv_c,
                                       const CartesianOrbitState &rv_d) {
    return CartesianOrbitState(Inertial2Synodic(rv_c.GetVec6(), rv_d.GetVec6()), rv_c.GetFrame());
  }

  Vec6 Inertial2Synodic(const Vec6 &rv_c, const Vec6 &rv_d) {
    Vec3 r_d = rv_d.head(3);
    Vec3 v_d = rv_d.tail(3);
    Vec3 r_c = rv_c.head(3);
    Vec3 v_c = rv_c.tail(3);

    Vec3 x = r_c.normalized();
    Vec3 y = (r_c.cross(v_c)).normalized();
    Vec3 z = y.cross(x);

    Mat3 R_inert2syn;
    R_inert2syn << x.transpose(), z.transpose(), y.transpose();

    Vec3 w = r_c.cross(v_c) / r_c.norm();
    Vec3 r_syn_d = R_inert2syn * (r_d - r_c);
    Vec3 v_syn_d = R_inert2syn * (v_d - v_c - w.cross(r_d - r_c));

    Vec6 rv_syn_d;
    rv_syn_d << r_syn_d, v_syn_d;
    return rv_syn_d;
  }

  CartesianOrbitState Synodic2Intertial(const CartesianOrbitState &rv_c,
                                        const CartesianOrbitState &rv_syn_d) {
    return CartesianOrbitState(Synodic2Intertial(rv_c.GetVec6(), rv_syn_d.GetVec6()),
                               rv_c.GetFrame());
  }

  Vec6 Synodic2Intertial(const Vec6 &rv_c, const Vec6 &rv_syn_d) {
    Vec3 r_c = rv_c.head(3);
    Vec3 v_c = rv_c.tail(3);

    // RTN basis vectors
    Vec3 x = r_c.normalized();
    Vec3 y = (r_c.cross(v_c)).normalized();
    Vec3 z = y.cross(x);

    Mat3 R_syn2inert;  // Rotation matrix from RTN to inertial
    R_syn2inert << x, z, y;

    Vec3 w = r_c.cross(v_c) / r_c.norm();
    Vec3 r_d = r_c + R_syn2inert * rv_syn_d.head(3);
    Vec3 v_d = v_c + R_syn2inert * rv_syn_d.tail(3) + w.cross(r_d - r_c);

    Vec6 rv_d;
    rv_d << r_d, v_d;
    return rv_d;
  }

  // From ClassicalOE
  // - To CartesianOrbitState
  CartesianOrbitState Classical2Cart(const ClassicalOE &coe, Real GM) {
    return CartesianOrbitState(Classical2Cart(coe.GetVec6(), GM), coe.GetFrame());
  }

  /// @brief Convert classical orbital elements to Cartesian state
  /// @param coe Classical orbital elements [km, –, rad, rad, rad, rad]
  /// @param GM Gravitational parameter [km^3/s^2]
  /// @return Cartesian state [km, km/s]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec6 Classical2Cart(const Vec6 &coe, Real GM) {
    auto [a, e, i, Omega, omega, M] = unpack(coe);

    Real E = Mean2EccAnomaly(M, e);  // Eccentric anomaly
    Real cosE = cos(E);
    Real sinE = sin(E);

    Real fac = sqrt((1.0 - e) * (1.0 + e));
    Real R = a * (1.0 - e * cosE);  // Distance
    Real V = sqrt(GM * a) / R;      // Velocity
    Vec3 r(a * (cosE - e), a * fac * sinE, 0.0);
    Vec3 v(-V * sinE, +V * fac * cosE, 0.0);

    Mat3 PQW = RotZ(-Omega) * RotX(-i) * RotZ(-omega);
    r = PQW * r;
    v = PQW * v;

    Vec6 rv;
    rv << r, v;
    return rv;
  }

  // - To QuasiNonsingOE
  QuasiNonsingOE Classical2QuasiNonsing(const ClassicalOE &coe, Real GM) {
    return QuasiNonsingOE(Classical2QuasiNonsing(coe.GetVec6(), GM), coe.GetFrame());
  }

  Vec6 Classical2QuasiNonsing(const Vec6 &coe, Real GM) {
    (void)GM;
    auto [a, e, i, Omega, w, M] = unpack(coe);

    Real u = w + M;
    Real ex = e * cos(w);
    Real ey = e * sin(w);

    return Vec6(a, u, ex, ey, i, Omega);
  }

  // - To EquinoctialOE
  EquinoctialOE Classical2Equinoctial(const ClassicalOE &coe, Real GM) {
    return EquinoctialOE(Classical2Equinoctial(coe.GetVec6(), GM), coe.GetFrame());
  }

  Vec6 Classical2Equinoctial(const Vec6 &coe, Real GM) {
    (void)GM;
    auto [a, e, i, Omega, w, M] = unpack(coe);

    Real f = Mean2TrueAnomaly(M, e);
    Real w_tilde = Omega + w;
    Real Psi = w_tilde + f;
    Real tq1 = e * cos(w_tilde);
    Real tq2 = e * sin(w_tilde);
    Real p1 = tan(i / 2) * cos(Omega);
    Real p2 = tan(i / 2) * sin(Omega);

    Psi = fmod(Psi.val(), 2 * PI);
    if (Psi > PI) {
      Psi = Psi - 2 * PI;
    }

    return Vec6(a, Psi, tq1, tq2, p1, p2);
  }

  // - To DelaunayOE
  DelaunayOE Classical2Delaunay(const ClassicalOE &coe, Real GM) {
    return DelaunayOE(Classical2Delaunay(coe.GetVec6(), GM), coe.GetFrame());
  }

  Vec6 Classical2Delaunay(const Vec6 &coe, Real GM) {
    auto [a, e, i, O, w, M] = unpack(coe);

    Real n = sqrt(GM / pow(a, 3));
    Real t = M / n;

    Real l = M;
    Real g = w;
    Real h = O - n * t;
    Real L = sqrt(GM * a);
    Real G = L * sqrt(1 - e * e);
    Real H = G * cos(i);

    return Vec6(l, g, h, L, G, H);
  }

  // From QuasiNonsingOE
  // - To ClassicalOE
  ClassicalOE QuasiNonsing2Classical(const QuasiNonsingOE &qnsoe, Real GM) {
    return ClassicalOE(QuasiNonsing2Classical(qnsoe.GetVec6(), GM), qnsoe.GetFrame());
  }

  Vec6 QuasiNonsing2Classical(const Vec6 &qnsoeVec, Real GM) {
    (void)GM;
    auto [a, u, ex, ey, i, Omega] = unpack(qnsoeVec);

    Real e = sqrt(ex * ex + ey * ey);
    Real w = atan2(ey, ex);
    Real M = u - w;

    Vec6 coe(a, e, i, Omega, w, M);
    return coe;
  }

  // From EquinoctialOE
  // - To ClassicalOE
  ClassicalOE Equinoctial2Classical(const EquinoctialOE &eqoe, Real GM) {
    return ClassicalOE(Equinoctial2Classical(eqoe.GetVec6(), GM), eqoe.GetFrame());
  }

  Vec6 Equinoctial2Classical(const Vec6 &equioe, Real GM) {
    (void)GM;
    auto [a, Psi, tq1, tq2, p1, p2] = unpack(equioe);

    Real Omega = atan2(p2, p1);
    Real i = 2 * atan2(p1, cos(Omega));

    Real wtilde = atan2(tq2, tq1);
    Real e = sqrt(tq1 * tq1 + tq2 * tq2);

    Real w = wtilde - Omega;
    Real f = Psi - wtilde;

    Real E = atan2(sin(f) * sqrt(1 - e * e), cos(f) + e);
    Real M = E - e * sin(E);

    w = std::fmod(w.val(), 2 * PI);
    M = std::fmod(M.val(), 2 * PI);
    if (std::abs(M.val() - 2 * PI) < std::numeric_limits<double>::epsilon()) {
      M = 0;
    }

    return Vec6(a, e, i, Omega, w, M);
  }

  // From DelaunayOE
  // - To ClassicalOE
  ClassicalOE Delaunay2Classical(const DelaunayOE &deloe, Real GM) {
    return ClassicalOE(Delaunay2Classical(deloe.GetVec6(), GM), deloe.GetFrame());
  }

  Vec6 Delaunay2Classical(const Vec6 &delaunay, Real GM) {
    auto [l, g, h, L, G, H] = unpack(delaunay);

    Real a = L * L / GM;
    Real M = l;

    Real n = sqrt(GM / pow(a, 3));
    Real t = M / n;

    Real e = sqrt(1 - pow(G / L, 2));
    Real i = safe_acos(H / G);
    Real O = h + n * t;
    Real w = g;

    return Vec6(a, e, i, O, w, M);
  }

  // From ClassicalOE QuasiNonsingOE
  // - To ClassicalOE
  ClassicalOE RelQuasiNonsing2Classical(const ClassicalOE &coe_c,
                                        const QuasiNonsingROE &RelQuasiNonsing) {
    return ClassicalOE(RelQuasiNonsing2Classical(coe_c.GetVec6(), RelQuasiNonsing.GetVec6()),
                       coe_c.GetFrame());
  }

  Vec6 RelQuasiNonsing2Classical(const Vec6 &coe_c, const Vec6 &RelQuasiNonsing) {
    auto [ac, ec, ic, Oc, wc, Mc] = unpack(coe_c);
    auto [ada, adl, adex, adey, adix, adiy] = unpack(RelQuasiNonsing);

    Real uc = Wrap2Pi(wc + Mc);
    Real exc = ec * cos(wc);
    Real eyc = ec * sin(wc);

    Real dO = (adiy / ac) / sin(ic);
    Real du = (adl / ac) - (dO * cos(ic));

    Real ad = ac + ada;
    Real exd = exc + (adex / ac);
    Real eyd = eyc + (adey / ac);
    Real id = ic + (wc / ac);
    Real Od = Oc + dO;
    Real ud = uc + du;

    Real wd = atan2(eyd, exd);
    Real ed = sqrt(exd * exd + eyd * eyd);
    Real Md = Wrap2Pi(ud - wd);

    return Vec6(ad, ed, id, Od, wd, Md);
  }

  VEC_IMP_VECTOR_REAL(Classical2Cart, 6);
  VEC_IMP_VECTOR_VECTOR(Inertial2Synodic, 6);
  VEC_IMP_VECTOR_VECTOR(Synodic2Intertial, 6);
  VEC_IMP_VECTOR_REAL(Cart2Classical, 6);
  VEC_IMP_VECTOR_REAL(Classical2QuasiNonsing, 6);
  VEC_IMP_VECTOR_REAL(Classical2Equinoctial, 6);
  VEC_IMP_VECTOR_REAL(Classical2Delaunay, 6);
  VEC_IMP_VECTOR_REAL(QuasiNonsing2Classical, 6);
  VEC_IMP_VECTOR_REAL(Equinoctial2Classical, 6);
  VEC_IMP_VECTOR_REAL(Delaunay2Classical, 6);
  VEC_IMP_VECTOR_VECTOR(RelQuasiNonsing2Classical, 6);

}  // namespace lupnt
