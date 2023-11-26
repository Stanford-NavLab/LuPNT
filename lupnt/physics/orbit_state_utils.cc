/**
 * @file OrbitStateUtils.cpp
 * @author Stanford NAV LAB
 * @brief Util functions for state conversions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "orbit_state_utils.h"

#include "lupnt/numerics/math_utils.h"

namespace lupnt {

CartesianOrbitState CoeToCart(const ClassicalOE &coe, double mu) {
  return CartesianOrbitState(CoeToCart(coe.GetVector(), mu),
                             coe.GetCoordSystem());
}

Vector6real CoeToCart(const Vector6real &coe, double mu) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  real p = a * (1.0 - pow(e, 2.0));
  real nu = MeanToTrueAnomaly(M, e);
  real pev = p / (1 + e * cos(nu));
  real mu_p = sqrt(mu / p);

  Vector3real r_PQW = {pev * cos(nu), pev * sin(nu), 0};
  Vector3real v_PQW = {-mu_p * sin(nu), mu_p * (e + cos(nu)), 0};

  Matrix3real rot{
      {cos(Omega) * cos(w) - sin(Omega) * sin(w) * cos(i),
       -cos(Omega) * sin(w) - sin(Omega) * cos(w) * cos(i),
       sin(Omega) * sin(i)},
      {sin(Omega) * cos(w) + cos(Omega) * sin(w) * cos(i),
       -sin(Omega) * sin(w) + cos(Omega) * cos(w) * cos(i),
       -cos(Omega) * sin(i)},
      {sin(w) * sin(i), cos(w) * sin(i), cos(i)},
  };

  Vector3real r = rot * r_PQW;
  Vector3real v = rot * v_PQW;

  Vector6real cart;
  cart << r, v;
  return cart;
}

ClassicalOE CartToCoe(const CartesianOrbitState &cart, double mu) {
  return ClassicalOE(CartToCoe(cart.GetVector(), mu), cart.GetCoordSystem());
}
Vector6real CartToCoe(const Vector6real &cart, double mu) {
  real a, e, p, i, Omega, w, nu, M;

  Vector3real r = cart.head(3);
  Vector3real v = cart.tail(3);

  real rnorm = r.squaredNorm();
  real vnorm = v.squaredNorm();

  Vector3real K{0., 0., 1.};

  Vector3real h = r.cross(v);
  Vector3real n = K.cross(h);

  real hnorm = h.squaredNorm();
  real nnorm = n.squaredNorm();

  Vector3real evec = ((pow(vnorm, 2.0) - mu / rnorm) * r - r.dot(v) * v) / mu;
  e = evec.norm();

  real xi = pow(vnorm, 2.0) / 2.0 - mu / rnorm;
  if (e.val() != 1.0) {
    a = -mu / (2.0 * xi);
    p = a * (1 - pow(e, 2));
  } else {
    a = INFINITY;
    p = pow(hnorm, 2.0) / mu;
  }

  i = acos(h(2) / hnorm);

  Omega = acos(n(0) / nnorm);
  if (n(1).val() < 0) {
    Omega = 2 * M_PI - Omega;
  }

  real ndote = n.dot(evec);
  w = acos(ndote / (nnorm * e));
  if (evec(2).val() < 0) {
    w = 2 * M_PI - w;
  }

  real edotr = evec.dot(r);
  if (edotr / e / rnorm >= 1.0) {
    nu = acos((1.0 - 1.0e-12) * abs(edotr / e / rnorm));
  } else if (edotr / e / rnorm <= -1.0) {
    nu = acos(-(1.0 - 1.0e-12) * abs(edotr / e / rnorm));
  } else {
    nu = acos(edotr / e / rnorm);  // true anomaly
  }
  if (edotr < 0) {
    nu = 2 * M_PI - nu;
  }

  M = TrueToMeanAnomaly(nu, e);

  return Vector6real{a, e, i, Omega, w, M};
}

/**
 * @brief Compute true anomaly from the eccentric anomaly
 *
 * @param E          Eccentric anomaly [rad]
 * @param e          Eccentricity
 * @return real  True anomaly [rad]
 */
real EccentricToTrueAnomaly(real E, real e) {
  return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e);
}

/**
 * @brief Compute mean anomaly from the eccentric anomaly
 *
 * @param E          Eccentric anomaly [rad]
 * @param e          Eccentricity
 * @return real  Mean anomaly [rad]
 */
real EccentricToMeanAnomaly(real E, real e) { return wrapToPi(E - e * sin(E)); }

/**
 * @brief Compute mean anomaly from the Eccentric anomaly
 * @param M          Mean anomaly [rad]
 * @param e          Eccentricity
 * @return real  Eccentric anomaly [rad]
 */
real MeanToEccentricAnomaly(real M, real e) {
  real MM = wrapToPi(M);

  // Initial estimate of E
  real E = MM;
  real Eest = E - (E - e * sin(E) - MM) / (1.0 - e * cos(E));

  double tol = 1e-9;
  int max_itr = 100;
  int itr = 0;

  while ((abs(Eest - E) >= tol) && (itr <= max_itr)) {
    E = Eest;
    Eest = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
    itr++;
  }
  E = Eest;

  return wrapToPi(E);
}

/**
 * @brief Compute mean anomaly from the Eccentric anomaly
 * @param M          Mean anomaly [rad]
 * @param e          Eccentricity
 * @return real  True anomaly [rad]
 */
real MeanToTrueAnomaly(real M, real e) {
  real E = MeanToEccentricAnomaly(M, e);
  real nu = EccentricToTrueAnomaly(E, e);
  return wrapToPi(nu);
}

/**
 * @brief Compute true anomaly from the Eccentric anomaly
 * @param nu         True anomaly [rad]
 * @param e          Eccentricity
 * @return real  Eccentric anomaly [rad]
 */
real TrueToEccentricAnomaly(real nu, real e) {
  real E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
  return E;
}

/**
 * @brief Compute mean anomaly from the True anomaly
 *
 * @param nu         True anomaly [rad]
 * @param e          Eccentricity
 * @return real  Mean anomaly [rad]
 */
real TrueToMeanAnomaly(real nu, real e) {
  real E = TrueToEccentricAnomaly(nu, e);
  real M = EccentricToMeanAnomaly(E, e);
  return M;
}

Vector6real QnsroeToCoe(const Vector6real &coe_c, const Vector6real &qnsroe) {
  auto [ac, ec, ic, Oc, wc, Mc] = unpack(coe_c);
  auto [ada, adl, adex, adey, adix, adiy] = unpack(qnsroe);

  real uc = wrapToPi(wc + Mc);
  real exc = ec * cos(wc);
  real eyc = ec * sin(wc);

  real dO = (adiy / ac) / sin(ic);
  real du = (adl / ac) - (dO * cos(ic));

  real ad = ac + ada;
  real exd = exc + (adex / ac);
  real eyd = eyc + (adey / ac);
  real id = ic + (wc / ac);
  real Od = Oc + dO;
  real ud = uc + du;

  real wd = atan2(eyd, exd);
  real ed = sqrt(exd * exd + eyd * eyd);
  real Md = wrapToPi(ud - wd);

  return Vector6real(ad, ed, id, Od, wd, Md);
}

ClassicalOE QnsroeToCoe(const ClassicalOE &coe_c,
                        const QuasiNonsingularROE &qnsroe) {
  return ClassicalOE(QnsroeToCoe(coe_c.GetVector(), qnsroe.GetVector()),
                     coe_c.GetCoordSystem());
}

Vector6real InertialToRtn(const Vector6real &rtn_c, const Vector6real &rtn_d) {
  Vector3real rInertial = rtn_d.head(3);
  Vector3real vInertial = rtn_d.tail(3);
  Vector3real rOrigin = rtn_c.head(3);
  Vector3real vOrigin = rtn_c.tail(3);

  Vector3real uR, uT, uN;  // RTN basis vectors
  uR = rOrigin.normalized();
  uN = (rOrigin.cross(vOrigin)).normalized();
  uT = uN.cross(uR);

  Matrix3real R_Inertial_Rtn;  // Rotation matrix from inertial to RTN
  R_Inertial_Rtn << uR.transpose(), uT.transpose(), uN.transpose();

  Vector3real w, rRtn, vRtn;
  w = rOrigin.cross(vOrigin) / rOrigin.squaredNorm();
  rRtn = R_Inertial_Rtn * (rInertial - rOrigin);
  vRtn = R_Inertial_Rtn * (vInertial - vOrigin - w.cross(rInertial - rOrigin));

  Vector6real rtnVec;
  rtnVec << rRtn, vRtn;
  return rtnVec;
}

CartesianOrbitState InertialToRtn(
    const CartesianOrbitState &rtn_c,
    const CartesianOrbitState &inertialOrbitState) {
  return CartesianOrbitState(
      InertialToRtn(rtn_c.GetVector(), inertialOrbitState.GetVector()),
      rtn_c.GetCoordSystem());
}

Vector6real CoeToRtn(const Vector6real &coe_c, const Vector6real &coe_d,
                     double mu) {
  Vector6real rtn_c = CoeToCart(coe_c, mu);
  Vector6real inertialVec = CoeToCart(coe_d, mu);
  Vector6real rtnVec = InertialToRtn(rtn_c, inertialVec);
  return rtnVec;
}

CartesianOrbitState CoeToRtn(const ClassicalOE &coe_c, const ClassicalOE &coe_d,
                             double mu) {
  Vector6real rtnVec = CoeToRtn(coe_c.GetVector(), coe_d.GetVector(), mu);
  return CartesianOrbitState(rtnVec);
}

QuasiNonsingularOE CoeToQnsoe(const ClassicalOE &coe) {
  Vector6real qnsoe = CoeToQnsoe(coe.GetVector());
  return QuasiNonsingularOE(qnsoe);
}

Vector6real CoeToQnsoe(const Vector6real &coe) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  real u = w + M;
  real ex = e * cos(w);
  real ey = e * sin(w);

  return Vector6real(a, u, ex, ey, i, Omega);
}

ClassicalOE QnsoeToCoe(const QuasiNonsingularOE &qnsoe) {
  return ClassicalOE(QnsoeToCoe(qnsoe.GetVector()), qnsoe.GetCoordSystem());
}
Vector6real QnsoeToCoe(const Vector6real &qnsoeVec) {
  auto [a, u, ex, ey, i, Omega] = unpack(qnsoeVec);

  real e = sqrt(ex * ex + ey * ey);
  real w = atan2(ey, ex);
  real M = u - w;

  Vector6real coe{a, e, i, Omega, w, M};
  return coe;
}

QuasiNonsingularROE QnsoeToQnsroe(const QuasiNonsingularOE &qnsoe_c,
                                  const QuasiNonsingularOE &qnsoe_d) {
  return QuasiNonsingularROE(
      QnsoeToQnsroe(qnsoe_c.GetVector(), qnsoe_d.GetVector()));
}

Vector6real QnsoeToQnsroe(const Vector6real &qnsoe_c,
                          const Vector6real &qnsoe_d) {
  auto [a_c, u_c, ex_c, ey_c, i_c, Omega_c] = unpack(qnsoe_c);
  auto [a_d, u_d, ex_d, ey_d, i_d, Omega_d] = unpack(qnsoe_d);

  real da, dl, dex, dey, dix, diy;
  da = (a_d - a_c) / a_c;
  dl = (u_d - u_c) + (Omega_d - Omega_c) * cos(i_c);
  dex = ex_d - ex_c;
  dey = ey_d - ey_c;
  dix = i_d - i_c;
  diy = (Omega_d - Omega_c) * sin(i_c);

  return a_c * Vector6real(da, dl, dex, dey, dix, diy);
}

QuasiNonsingularROE CoeToQnsroe(const ClassicalOE &coe_c,
                                const ClassicalOE &coe_d) {
  return QuasiNonsingularROE(
      QnsoeToQnsroe(CoeToQnsoe(coe_c), CoeToQnsoe(coe_d)));
}

EquinoctialOE CoeToEqoe(const ClassicalOE &coe) {
  return EquinoctialOE(CoeToEqoe(coe.GetVector()), coe.GetCoordSystem());
}
ClassicalOE EqoeToCoe(const EquinoctialOE &eqoe) {
  return ClassicalOE(EqoeToCoe(eqoe.GetVector()), eqoe.GetCoordSystem());
}

Vector6real EqoeToCoe(const Vector6real &equioe) {
  auto [a, Psi, tq1, tq2, p1, p2] = unpack(equioe);

  real Omega = atan2(p2, p1);
  real i = 2 * atan2(p1, cos(Omega));

  real wtilde = atan2(tq2, tq1);
  real e = sqrt(tq1 * tq1 + tq2 * tq2);

  real w = wtilde - Omega;
  real f = Psi - wtilde;

  real E = atan2(sin(f) * sqrt(1 - e * e), cos(f) + e);
  real M = E - e * sin(E);

  w = std::fmod(w.val(), 2 * M_PI);
  M = std::fmod(M.val(), 2 * M_PI);
  if (std::abs(M.val() - 2 * M_PI) < std::numeric_limits<double>::epsilon()) {
    M = 0;
  }

  return Vector6real(a, e, i, Omega, w, M);
}

Vector6real MeanToOsc(const Vector6real &coe_m, double J2) {
  Vector6real coe_o;

  if (J2 > 0) {
    Vector6real meanEquioe = CoeToEqoe(coe_m);
    Vector6real oscEquioe;  // = MeanOscClosedEqui(meanEquioe, J2);
    coe_o = EqoeToCoe(oscEquioe);
  } else {
    coe_o = coe_m;
  }

  return coe_o;
}

ClassicalOE MeanToOsc(const ClassicalOE &coe_m, double J2) {
  return ClassicalOE(MeanToOsc(coe_m.GetVector(), J2), coe_m.GetCoordSystem());
}

Vector6real osc2mean_NRiterator(const Vector6real &osc_equi_elem, double tol) {
  Vector6real mean_equi_elem = osc_equi_elem;
  double R = 1.0;
  int niter = 0;

  while (std::abs(R) > tol) {
    niter++;
    Vector6real osc_loop;
    // std::tie(std::ignore, osc_loop, std::ignore) =
    // transformationmatrix_osc2mean_equinoctial(mean_equi_elem);
    Vector6real delta = osc_equi_elem - osc_loop;
    // R = delta.norm_inf(); // Assuming the library provides an infinity norm
    // function
    mean_equi_elem = mean_equi_elem + delta;

    if (niter > 100) {
      std::cout << "Osc2Mean iterations > 100" << std::endl;
      break;
    }
  }

  return mean_equi_elem;
}

Vector6real OscToMean(const Vector6real &coe_o, double J2) {
  Vector6real coe_m;
  double tol = 1e-8;

  if (J2 > 0) {
    Vector6real eqoe_o = CoeToEqoe(coe_o);
    Vector6real eqoe_m = osc2mean_NRiterator(eqoe_o, tol);
    coe_m = EqoeToCoe(eqoe_m);
  } else {
    coe_m = coe_o;
  }

  return coe_m;
}

ClassicalOE OscToMean(const ClassicalOE &coe_o, double J2) {
  return ClassicalOE(OscToMean(coe_o.GetVector(), J2), coe_o.GetCoordSystem());
}

Vector6real CoeToEqoe(const Vector6real &coe) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  real f = MeanToTrueAnomaly(M, e);
  real w_tilde = Omega + w;
  real Psi = w_tilde + f;
  real tq1 = e * cos(w_tilde);
  real tq2 = e * sin(w_tilde);
  real p1 = tan(i / 2) * cos(Omega);
  real p2 = tan(i / 2) * sin(Omega);

  Psi = fmod(Psi.val(), 2 * M_PI);
  if (Psi > M_PI) {
    Psi = Psi - 2 * M_PI;
  }

  return Vector6real(a, Psi, tq1, tq2, p1, p2);
}

Vector6real CartToQnsoe(const Vector6real &cart, double mu) {
  Vector6real coe = CartToCoe(cart, mu);
  Vector6real qnsoe = CoeToQnsoe(coe);
  return qnsoe;
}
QuasiNonsingularOE CartToQnsoe(const CartesianOrbitState &cart, double mu) {
  return QuasiNonsingularOE(CartToQnsoe(cart.GetVector(), mu),
                            cart.GetCoordSystem());
}

Vector6real DelaunayToCoe(const Vector6real &delaunay, double mu, double n,
                          double t) {
  auto [l, g, h, L, G, H] = unpack(delaunay);

  real a = L * L / mu;
  real e = sqrt(1 - pow(G / L, 2));
  real i = acos(H / G);
  real O = h + n * t;
  real w = g;
  real M = l;

  return Vector6real(a, e, i, O, w, M);
}

Vector6real CoeToDelaunay(const Vector6real &coe, double mu, double n,
                          double t) {
  auto [a, e, i, O, w, M] = unpack(coe);

  real l = M;
  real g = w;
  real h = O - n * t;
  real L = sqrt(mu * a);
  real G = L * sqrt(1 - e * e);
  real H = G * cos(i);

  return Vector6real(l, g, h, L, G, H);
}

}  // namespace lupnt