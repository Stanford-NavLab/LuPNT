#include "conversions.h"

#include <cmath>

#include "anomaly.h"
#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "orbit_states.h"

namespace lupnt {
// From CartesianOrbitState
// - To ClassicalOE
ClassicalOE Cart2Classical(const CartesianOrbitState &rv, Real mu) {
  return ClassicalOE(Cart2Classical(rv.GetVec(), mu), rv.GetCoordSystem());
}

Vec6 Cart2Classical(const Vec6 &rv, Real mu) {
  Vec3 r = rv.head(3);
  Vec3 v = rv.tail(3);

  Vec3 k(0, 0, 1);
  Vec3 h = r.cross(v);
  Vec3 n = k.cross(h);

  Vec3 evec = v.cross(h) / mu - r / r.norm();

  Real e = evec.norm();
  Real i = safe_acos(h(2) / h.norm());
  Real a = 1.0 / (2.0 / r.norm() - pow(v.norm(), 2.0) / mu);

  Real nu =
      (r.dot(v) >= 0) ? angleBetweenVecs(evec, r) : -angleBetweenVecs(evec, r);

  Real Omega =
      (n(1) >= 0) ? safe_acos(n(0) / n.norm()) : -safe_acos(n(0) / n.norm());

  Real w =
      (evec(2) >= 0) ? angleBetweenVecs(n, evec) : -angleBetweenVecs(n, evec);

  Real M = True2MeanAnomaly(nu, e);

  return Vec6{a, e, i, Omega, w, M};
}

// - To CartesianOrbitState (relative)
CartesianOrbitState Inertial2Rtn(const CartesianOrbitState &rv_c,
                                 const CartesianOrbitState &rv_d) {
  return CartesianOrbitState(Inertial2Rtn(rv_c.GetVec(), rv_d.GetVec()),
                             rv_c.GetCoordSystem());
}

Vec6 Inertial2Rtn(const Vec6 &rv_c, const Vec6 &rv_d) {
  Vec3 r_d = rv_d.head(3);
  Vec3 v_d = rv_d.tail(3);
  Vec3 r_c = rv_c.head(3);
  Vec3 v_c = rv_c.tail(3);

  // RTN basis vectors
  Vec3 uR = r_c.normalized();
  Vec3 uN = (r_c.cross(v_c)).normalized();
  Vec3 uT = uN.cross(uR);

  Mat3 Rot_inert_rtn;  // Rotation matrix from inertial to RTN
  Rot_inert_rtn << uR.transpose(), uT.transpose(), uN.transpose();

  Vec3 w = r_c.cross(v_c) / r_c.norm();
  Vec3 r_rtn_d = Rot_inert_rtn * (r_d - r_c);
  Vec3 v_rtn_d = Rot_inert_rtn * (v_d - v_c - w.cross(r_d - r_c));

  Vec6 rv_rtn_d;
  rv_rtn_d << r_rtn_d, v_rtn_d;
  return rv_rtn_d;
}

CartesianOrbitState Rtn2Inertial(const CartesianOrbitState &rv_c,
                                 const CartesianOrbitState &rv_rtn_d) {
  return CartesianOrbitState(Rtn2Inertial(rv_c.GetVec(), rv_rtn_d.GetVec()),
                             rv_c.GetCoordSystem());
}

Vec6 Rtn2Inertial(const Vec6 &rv_c, const Vec6 &rv_rtn_d) {
  Vec3 r_c = rv_c.head(3);
  Vec3 v_c = rv_c.tail(3);

  // RTN basis vectors
  Vec3 uR = r_c.normalized();
  Vec3 uN = (r_c.cross(v_c)).normalized();
  Vec3 uT = uN.cross(uR);

  Mat3 Rot_rtn_inert;  // Rotation matrix from RTN to inertial
  Rot_rtn_inert << uR, uT, uN;

  Vec3 w = r_c.cross(v_c) / r_c.norm();
  Vec3 r_d = r_c + Rot_rtn_inert * rv_rtn_d.head(3);
  Vec3 v_d = v_c + Rot_rtn_inert * rv_rtn_d.tail(3) + w.cross(r_d - r_c);

  Vec6 rv_d;
  rv_d << r_d, v_d;
  return rv_d;
}

// From ClassicalOE
// - To CartesianOrbitState
CartesianOrbitState Classical2Cart(const ClassicalOE &coe, Real mu) {
  return CartesianOrbitState(Classical2Cart(coe.GetVec(), mu),
                             coe.GetCoordSystem());
}

Vec6 Classical2Cart(const Vec6 &coe, Real mu) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  Real p = a * (1.0 - pow(e, 2.0));
  Real nu = Mean2TrueAnomaly(M, e);
  Real pev = p / (1.0 + e * cos(nu));
  Real mu_p = sqrt(mu / p);

  Vec3 r_PQW{pev * cos(nu), pev * sin(nu), 0.0};
  Vec3 v_PQW{-mu_p * sin(nu), mu_p * (e + cos(nu)), 0.0};

  // rot = RotZ(-Omega) * RotX(-i) * RotZ(-w)
  Mat3 rot{
      {cos(Omega) * cos(w) - sin(Omega) * sin(w) * cos(i),
       -cos(Omega) * sin(w) - sin(Omega) * cos(w) * cos(i),
       sin(Omega) * sin(i)},
      {sin(Omega) * cos(w) + cos(Omega) * sin(w) * cos(i),
       -sin(Omega) * sin(w) + cos(Omega) * cos(w) * cos(i),
       -cos(Omega) * sin(i)},
      {sin(w) * sin(i), cos(w) * sin(i), cos(i)},
  };

  Vec3 r = rot * r_PQW;
  Vec3 v = rot * v_PQW;

  Vec6 rv;
  rv << r, v;
  return rv;
}

// - To QuasiNonsingOE
QuasiNonsingOE Classical2QuasiNonsing(const ClassicalOE &coe, Real mu) {
  return QuasiNonsingOE(Classical2QuasiNonsing(coe.GetVec(), mu),
                        coe.GetCoordSystem());
}

Vec6 Classical2QuasiNonsing(const Vec6 &coe, Real mu) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  Real u = w + M;
  Real ex = e * cos(w);
  Real ey = e * sin(w);

  return Vec6(a, u, ex, ey, i, Omega);
}

// - To EquinoctialOE
EquinoctialOE Classical2Equinoctial(const ClassicalOE &coe, Real mu) {
  return EquinoctialOE(Classical2Equinoctial(coe.GetVec(), mu),
                       coe.GetCoordSystem());
}

Vec6 Classical2Equinoctial(const Vec6 &coe, Real mu) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  Real f = Mean2TrueAnomaly(M, e);
  Real w_tilde = Omega + w;
  Real Psi = w_tilde + f;
  Real tq1 = e * cos(w_tilde);
  Real tq2 = e * sin(w_tilde);
  Real p1 = tan(i / 2) * cos(Omega);
  Real p2 = tan(i / 2) * sin(Omega);

  Psi = fmod(Psi.val(), 2 * M_PI);
  if (Psi > M_PI) {
    Psi = Psi - 2 * M_PI;
  }

  return Vec6(a, Psi, tq1, tq2, p1, p2);
}

// - To DelaunayOE
DelaunayOE Classical2Delaunay(const ClassicalOE &coe, Real mu) {
  return DelaunayOE(Classical2Delaunay(coe.GetVec(), mu), coe.GetCoordSystem());
}

Vec6 Classical2Delaunay(const Vec6 &coe, Real mu) {
  auto [a, e, i, O, w, M] = unpack(coe);

  Real n = sqrt(mu / pow(a, 3));
  Real t = M / n;

  Real l = M;
  Real g = w;
  Real h = O - n * t;
  Real L = sqrt(mu * a);
  Real G = L * sqrt(1 - e * e);
  Real H = G * cos(i);

  return Vec6(l, g, h, L, G, H);
}

// From QuasiNonsingOE
// - To ClassicalOE
ClassicalOE QuasiNonsing2Classical(const QuasiNonsingOE &qnsoe, Real mu) {
  return ClassicalOE(QuasiNonsing2Classical(qnsoe.GetVec(), mu),
                     qnsoe.GetCoordSystem());
}

Vec6 QuasiNonsing2Classical(const Vec6 &qnsoeVec, Real mu) {
  auto [a, u, ex, ey, i, Omega] = unpack(qnsoeVec);

  Real e = sqrt(ex * ex + ey * ey);
  Real w = atan2(ey, ex);
  Real M = u - w;

  Vec6 coe{a, e, i, Omega, w, M};
  return coe;
}

// From EquinoctialOE
// - To ClassicalOE
ClassicalOE Equinoctial2Classical(const EquinoctialOE &eqoe, Real mu) {
  return ClassicalOE(Equinoctial2Classical(eqoe.GetVec(), mu),
                     eqoe.GetCoordSystem());
}

Vec6 Equinoctial2Classical(const Vec6 &equioe, Real mu) {
  auto [a, Psi, tq1, tq2, p1, p2] = unpack(equioe);

  Real Omega = atan2(p2, p1);
  Real i = 2 * atan2(p1, cos(Omega));

  Real wtilde = atan2(tq2, tq1);
  Real e = sqrt(tq1 * tq1 + tq2 * tq2);

  Real w = wtilde - Omega;
  Real f = Psi - wtilde;

  Real E = atan2(sin(f) * sqrt(1 - e * e), cos(f) + e);
  Real M = E - e * sin(E);

  w = std::fmod(w.val(), 2 * M_PI);
  M = std::fmod(M.val(), 2 * M_PI);
  if (std::abs(M.val() - 2 * M_PI) < std::numeric_limits<double>::epsilon()) {
    M = 0;
  }

  return Vec6(a, e, i, Omega, w, M);
}

// From DelaunayOE
// - To ClassicalOE
ClassicalOE Delaunay2Classical(const DelaunayOE &deloe, Real mu) {
  return ClassicalOE(Delaunay2Classical(deloe.GetVec(), mu),
                     deloe.GetCoordSystem());
}

Vec6 Delaunay2Classical(const Vec6 &delaunay, Real mu) {
  auto [l, g, h, L, G, H] = unpack(delaunay);

  Real a = L * L / mu;
  Real M = l;

  Real n = sqrt(mu / pow(a, 3));
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
  return ClassicalOE(
      RelQuasiNonsing2Classical(coe_c.GetVec(), RelQuasiNonsing.GetVec()),
      coe_c.GetCoordSystem());
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
VEC_IMP_VECTOR_VECTOR(Inertial2Rtn, 6);
VEC_IMP_VECTOR_VECTOR(Rtn2Inertial, 6);
VEC_IMP_VECTOR_REAL(Cart2Classical, 6);
VEC_IMP_VECTOR_REAL(Classical2QuasiNonsing, 6);
VEC_IMP_VECTOR_REAL(Classical2Equinoctial, 6);
VEC_IMP_VECTOR_REAL(Classical2Delaunay, 6);
VEC_IMP_VECTOR_REAL(QuasiNonsing2Classical, 6);
VEC_IMP_VECTOR_REAL(Equinoctial2Classical, 6);
VEC_IMP_VECTOR_REAL(Delaunay2Classical, 6);
VEC_IMP_VECTOR_VECTOR(RelQuasiNonsing2Classical, 6);

}  // namespace lupnt