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

/**
 * @brief Convert classical orbital elements to Cartesian
 *
 * @param coe
 * @return Vector6real
 * @ref Vallado "Fudamentals of Astrodynamics and Applications " p146 (ELORB)
 */
CartesianOrbitState CoeToCart(const ClassicalOE &coe, double mu) {
  Vector6real coeVec = coe.GetVector();
  Vector6real cartVec = CoeToCart(coeVec, mu);
  CartesianOrbitState cartOrbitState(cartVec);
  cartOrbitState.SetCoordSystem(coe.GetCoordSystem());
  return cartOrbitState;
}

Vector6real CoeToCart(const Vector6real &coeVec, double mu) {
  real a = coeVec(0);
  real e = coeVec(1);
  real i = coeVec(2);
  real Omega = coeVec(3);
  real w = coeVec(4);
  real M = coeVec(5);

  real p = a * (1.0 - pow(e, 2.0));
  real nu = MeanAnomToTrueAnom(M, e);
  real pev = p / (1 + e * cos(nu));
  real mu_p = sqrt(mu / p);

  Vector3real r_PQW = {pev * cos(nu), pev * sin(nu), 0};
  Vector3real v_PQW = {-mu_p * sin(nu), mu_p * (e + cos(nu)), 0};

  Matrix3real rot;
  rot.row(0) << cos(Omega) * cos(w) - sin(Omega) * sin(w) * cos(i),
      -cos(Omega) * sin(w) - sin(Omega) * cos(w) * cos(i), sin(Omega) * sin(i);
  rot.row(1) << sin(Omega) * cos(w) + cos(Omega) * sin(w) * cos(i),
      -sin(Omega) * sin(w) + cos(Omega) * cos(w) * cos(i), -cos(Omega) * sin(i);
  rot.row(2) << sin(w) * sin(i), cos(w) * sin(i), cos(i);

  Vector3real r = rot * r_PQW;
  Vector3real v = rot * v_PQW;

  Vector6real cartVec;
  cartVec << r, v;
  return cartVec;
}

/**
 * @brief Convert Cartesian to classical orbital elements
 *
 * @param x_cart
 * @return Vector6real
 * @ref Vallado "Fudamentals of Astrodynamics and Applications " p146 (ELORB)
 */
ClassicalOE CartToCoe(const CartesianOrbitState &cartOrbitState, double mu) {
  Vector6real cartVec = cartOrbitState.GetVector();
  Vector6real coeVec = CartToCoe(cartVec, mu);
  ClassicalOE coe(coeVec);
  return coe;
}
Vector6real CartToCoe(const Vector6real &cartVec, double mu) {
  real a, e, p, i, Omega, w, nu, M;

  Vector3real r = cartVec.head(3);
  Vector3real v = cartVec.tail(3);

  real rnorm = r.squaredNorm();
  real vnorm = v.squaredNorm();

  Vector3real K;
  K << 0, 0, 1.0;

  Vector3real h;
  h = r.cross(v);
  Vector3real n;
  n = K.cross(h);

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

  M = TrueAnomToMeanAnom(nu, e);

  Vector6real coeVec;
  coeVec << a, e, i, Omega, w, M;
  return coeVec;
}

/**
 * @brief Compute true anomaly from the eccentric anomaly
 *
 * @param E          Eccentric anomaly [rad]
 * @param e          Eccentricity
 * @return real  True anomaly [rad]
 */

real EccentricAnomToTrueAnom(real E, real e) {
  return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e);
  // real beta = e / (1.0 + sqrt(1.0 - pow(e, 2)));
  // return E + 2.0 * atan2(beta * sin(E), 1.0 + beta * cos(E));
}

/**
 * @brief Compute mean anomaly from the eccentric anomaly
 *
 * @param E          Eccentric anomaly [rad]
 * @param e          Eccentricity
 * @return real  Mean anomaly [rad]
 */

real EccentricAnomToMeanAnom(real E, real e) {
  return wrapToPi(E - e * sin(E));
}

/**
 * @brief Compute mean anomaly from the Eccentric anomaly
 * @param M          Mean anomaly [rad]
 * @param e          Eccentricity
 * @return real  Eccentric anomaly [rad]
 */

real MeanAnomToEccentricAnom(real M, real e) {
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

real MeanAnomToTrueAnom(real M, real e) {
  real E = MeanAnomToEccentricAnom(M, e);
  real nu = EccentricAnomToTrueAnom(E, e);
  return wrapToPi(nu);
}

/**
 * @brief Compute true anomaly from the Eccentric anomaly
 * @param nu         True anomaly [rad]
 * @param e          Eccentricity
 * @return real  Eccentric anomaly [rad]
 */

real TrueAnomToEccentricAnom(real nu, real e) {
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

real TrueAnomToMeanAnom(real nu, real e) {
  real E = TrueAnomToEccentricAnom(nu, e);
  real M = EccentricAnomToMeanAnom(E, e);
  return M;
}

Vector6real RoeToCoe(const Vector6real &coe_c, const Vector6real &roe) {
  real ac = coe_c[0];
  real ec = coe_c[1];
  real ic = coe_c[2];
  real Oc = coe_c[3];
  real uc = wrapToPi(coe_c[4] + coe_c[5]);
  real exc = ec * cos(coe_c[4]);
  real eyc = ec * sin(coe_c[4]);

  real dO = roe[5] / (ac * sin(ic));
  real du = (roe[1] / ac) - (dO * cos(ic));

  real ad = ac + roe[0];
  real exd = exc + (roe[2] / ac);
  real eyd = eyc + (roe[3] / ac);
  real id = ic + (roe[4] / ac);
  real Od = Oc + dO;
  real ud = uc + du;

  real wd = atan2(eyd, exd);
  real ed = sqrt(exd * exd + eyd * eyd);
  real Md = wrapToPi(ud - wd);

  Vector6real coe_d;
  coe_d << ad, ed, id, Od, wd, Md;

  return coe_d;
}

ClassicalOE RoeToCoe(const ClassicalOE &coe_c, const QuasiNonsingularROE &roe) {
  Vector6real coe_d = RoeToCoe(coe_c.GetVector(), roe.GetVector());
  return ClassicalOE(coe_d);
}

Vector6real InertialToRtn(const Vector6real &rtnOrigin,
                          const Vector6real &inertialVec) {
  Vector3real rInertial = inertialVec.head(3);
  Vector3real vInertial = inertialVec.tail(3);
  Vector3real rOrigin = rtnOrigin.head(3);
  Vector3real vOrigin = rtnOrigin.tail(3);

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
    const CartesianOrbitState &rtnOrigin,
    const CartesianOrbitState &inertialOrbitState) {
  Vector6real rtnOriginVec = rtnOrigin.GetVector();
  Vector6real inertialVec = inertialOrbitState.GetVector();
  Vector6real rtnVec = InertialToRtn(rtnOriginVec, inertialVec);
  CartesianOrbitState rtnOrbitState(rtnVec);
  return rtnOrbitState;
}

Vector6real CoeToRtn(const Vector6real &coe_c, const Vector6real &coe_d,
                     double mu) {
  Vector6real rtnOrigin = CoeToCart(coe_c, mu);
  Vector6real inertialVec = CoeToCart(coe_d, mu);
  Vector6real rtnVec = InertialToRtn(rtnOrigin, inertialVec);
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

Vector6real CoeToQnsoe(const Vector6real &coeVec) {
  real a = coeVec(0);
  real e = coeVec(1);
  real i = coeVec(2);
  real Omega = coeVec(3);
  real w = coeVec(4);
  real M = coeVec(5);

  real u = w + M;
  real ex = e * cos(w);
  real ey = e * sin(w);

  Vector6real qnsoeVec;
  qnsoeVec << a, u, ex, ey, i, Omega;
  return qnsoeVec;
}
ClassicalOE QnsoeToCoe(const QuasiNonsingularOE &qnsoe) {
  Vector6real coeVec = QnsoeToCoe(qnsoe.GetVector());
  ClassicalOE coe(coeVec);
  return coe;
}
Vector6real QnsoeToCoe(const Vector6real &qnsoeVec) {
  real a = qnsoeVec(0);
  real u = qnsoeVec(1);
  real ex = qnsoeVec(2);
  real ey = qnsoeVec(3);
  real i = qnsoeVec(4);
  real Omega = qnsoeVec(5);

  real e = sqrt(ex * ex + ey * ey);
  real w = atan2(ey, ex);
  real M = u - w;

  Vector6real coeVec;
  coeVec << a, e, i, Omega, w, M;
  return coeVec;
}

QuasiNonsingularROE QnsoeToQnsroe(const QuasiNonsingularOE &qnsoe_c,
                                  const QuasiNonsingularOE &qnsoe_d) {
  real a_c = qnsoe_c.a();
  real u_c = qnsoe_c.u();
  real ex_c = qnsoe_c.ex();
  real ey_c = qnsoe_c.ey();
  real i_c = qnsoe_c.i();
  real Omega_c = qnsoe_c.Omega();

  real a_d = qnsoe_d.a();
  real u_d = qnsoe_d.u();
  real ex_d = qnsoe_d.ex();
  real ey_d = qnsoe_d.ey();
  real i_d = qnsoe_d.i();
  real Omega_d = qnsoe_d.Omega();

  real da = (a_d - a_c) / a_c;
  real dl = (u_d - u_c) + (Omega_d - Omega_c) * cos(i_c);
  real dex = ex_d - ex_c;
  real dey = ey_d - ey_c;
  real dix = i_d - i_c;
  real diy = (Omega_d - Omega_c) * sin(i_c);

  Vector6real qnsroe;
  qnsroe << da, dl, dex, dey, dix, diy;
  return QuasiNonsingularROE(qnsroe);
}

Vector6real QnsoeToQnsroe(const Vector6real &qnsoe_c,
                          const Vector6real &qnsoe_d) {
  QuasiNonsingularROE qnsroe =
      QnsoeToQnsroe(QuasiNonsingularOE(qnsoe_c), QuasiNonsingularOE(qnsoe_d));
  return qnsroe.GetVector();
}

QuasiNonsingularROE CoeToQnsroe(const ClassicalOE &coe_c,
                                const ClassicalOE &coe_d) {
  QuasiNonsingularOE qnsoe_c = CoeToQnsoe(coe_c);
  QuasiNonsingularOE qnsoe_d = CoeToQnsoe(coe_d);
  QuasiNonsingularROE qnsroe = QnsoeToQnsroe(qnsoe_c, qnsoe_d);
  return qnsroe;
}

Vector6real EquioeToCoe(const Vector6real &equioe) {
  real a = equioe(0);
  real Psi = equioe(1);
  real tq1 = equioe(2);
  real tq2 = equioe(3);
  real p1 = equioe(4);
  real p2 = equioe(5);

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

  Vector6real coe;
  coe << a, e, i, Omega, w, M;

  return coe;
}

Vector6real MeanToOsculating(const Vector6real &meanCoe, double J2) {
  Vector6real oscCoe;

  if (J2 > 0) {
    Vector6real meanEquioe = CoeToEquioe(meanCoe);
    Vector6real oscEquioe;  // = MeanOscClosedEqui(meanEquioe, J2);
    oscCoe = EquioeToCoe(oscEquioe);
  } else {
    oscCoe = meanCoe;
  }

  return oscCoe;
}

ClassicalOE MeanToOsculating(const ClassicalOE &meanCoe, double J2) {
  Vector6real meanCoeVec = meanCoe.GetVector();
  Vector6real oscCoeVec = MeanToOsculating(meanCoeVec, J2);
  ClassicalOE oscCoe(oscCoeVec);
  return oscCoe;
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

Vector6real OsculatingToMean(const Vector6real &oscCoe, double J2) {
  Vector6real meanCoe;
  double tol = 1e-8;

  if (J2 > 0) {
    Vector6real oscEquioe = CoeToEquioe(oscCoe);
    Vector6real meanEquioe = osc2mean_NRiterator(oscEquioe, tol);
    meanCoe = EquioeToCoe(meanEquioe);
  } else {
    meanCoe = oscCoe;
  }

  return meanCoe;
}

ClassicalOE OsculatingToMean(const ClassicalOE &oscCoe, double J2) {
  Vector6real oscCoeVec = oscCoe.GetVector();
  Vector6real meanCoeVec = OsculatingToMean(oscCoeVec, J2);
  ClassicalOE meanCoe(meanCoeVec);
  return meanCoe;
}

Vector6real CoeToEquioe(const Vector6real &coe) {
  real a = coe(0);
  real e = coe(1);
  real i = coe(2);
  real Omega = coe(3);
  real w = coe(4);
  real M = coe(5);

  real f = MeanAnomToTrueAnom(M, e);
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

  Vector6real equioe;
  equioe << a, Psi, tq1, tq2, p1, p2;

  return equioe;
}

Vector6real CartToQnsoe(const Vector6real &cart, double mu) {
  Vector6real coe = CartToCoe(cart, mu);
  Vector6real qnsoe = CoeToQnsoe(coe);
  return qnsoe;
}
QuasiNonsingularOE CartToQnsoe(const CartesianOrbitState &cart, double mu) {
  Vector6real cartVec = cart.GetVector();
  Vector6real qnsoeVec = CartToQnsoe(cartVec, mu);
  QuasiNonsingularOE qnsoe(qnsoeVec);
  return qnsoe;
}

QuasiNonsingularOE OscQnsoeToMeanQnsoe(const QuasiNonsingularOE &oscQnsoe,
                                       double J2) {
  ClassicalOE oscCoe = QnsoeToCoe(oscQnsoe);
  ClassicalOE meanCoe = OsculatingToMean(oscCoe, J2);
  QuasiNonsingularOE meanQnsoe = CoeToQnsoe(meanCoe);
  return meanQnsoe;
}

Vector6real DelaunayToCoe(const Vector6real &delaunay, double mu, double n,
                          double t) {
  real l = delaunay(0);
  real g = delaunay(1);
  real h = delaunay(2);
  real L = delaunay(3);
  real G = delaunay(4);
  real H = delaunay(5);

  real a = L * L / mu;
  real e = sqrt(1 - pow(G / L, 2));
  real i = acos(H / G);
  real O = h + n * t;
  real w = g;
  real M = l;

  Vector6real coe;
  coe << a, e, i, O, w, M;
  return coe;
}

Vector6real CoeToDelaunay(const Vector6real &coe, double mu, double n,
                          double t) {
  real a = coe(0);
  real e = coe(1);
  real i = coe(2);
  real O = coe(3);
  real w = coe(4);
  real M = coe(5);

  real l = M;
  real g = w;
  real h = O - n * t;
  real L = sqrt(mu * a);
  real G = L * sqrt(1 - e * e);
  real H = G * cos(i);

  Vector6real delaunay;
  delaunay << l, g, h, L, G, H;
  return delaunay;
}

}  // namespace lupnt