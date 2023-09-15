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

#include "lupnt/physics/OrbitStateUtils.h"

#include <Eigen/Dense>

#include "lupnt/core/Constants.h"
#include "lupnt/numerics/MathUtils.h"

namespace LPT {

/* *********************************************************************************************
 */
/* **********************************  OrbitStateUtils
 * ********************************************** */
/* *********************************************************************************************
 */

/**
 * @brief Convert classical orbital elements to Cartesian
 *
 * @param coe
 * @return ad::Vector6real
 * @ref Vallado "Fudamentals of Astrodynamics and Applications " p146 (ELORB)
 */
CartesianOrbitState CoeToCart(const ClassicalOE coe, double mu) {
  ad::Vector6real coeVec = coe.GetVector();
  ad::Vector6real cartVec = CoeToCart(coeVec, mu);
  CartesianOrbitState cartOrbitState(cartVec);
  cartOrbitState.SetCoordSystem(coe.GetCoordSystem());
  return cartOrbitState;
}
template <typename T>
Eigen::Matrix<T, 6, 1> CoeToCart(const Eigen::Matrix<T, 6, 1> &coeVec,
                                 double mu) {
  T a = coeVec(0);
  T e = coeVec(1);
  T i = coeVec(2);
  T Omega = coeVec(3);
  T w = coeVec(4);
  T M = coeVec(5);

  T p = a * (1.0 - pow(e, 2.0));
  T nu = MeanAnomToTrueAnom(M, e);
  T pev = p / (1 + e * cos(nu));
  T mu_p = sqrt(mu / p);

  Eigen::Matrix<T, 3, 1> r_PQW;
  r_PQW << pev * cos(nu), pev * sin(nu), 0;
  Eigen::Matrix<T, 3, 1> v_PQW;
  v_PQW << -mu_p * sin(nu), mu_p * (e + cos(nu)), 0;

  ad::MatrixXreal rot(3, 3);
  rot.row(0) << cos(Omega) * cos(w) - sin(Omega) * sin(w) * cos(i),
      -cos(Omega) * sin(w) - sin(Omega) * cos(w) * cos(i), sin(Omega) * sin(i);
  rot.row(1) << sin(Omega) * cos(w) + cos(Omega) * sin(w) * cos(i),
      -sin(Omega) * sin(w) + cos(Omega) * cos(w) * cos(i), -cos(Omega) * sin(i);
  rot.row(2) << sin(w) * sin(i), cos(w) * sin(i), cos(i);

  Eigen::Matrix<T, 3, 1> r = rot * r_PQW;
  Eigen::Matrix<T, 3, 1> v = rot * v_PQW;

  Eigen::Matrix<T, 6, 1> cartVec;
  cartVec << r, v;
  return cartVec;
}

/**
 * @brief Convert Cartesian to classical orbital elements
 *
 * @param x_cart
 * @return ad::Vector6real
 * @ref Vallado "Fudamentals of Astrodynamics and Applications " p146 (ELORB)
 */
ClassicalOE CartToCoe(const CartesianOrbitState cartOrbitState, double mu) {
  ad::Vector6real cartVec = cartOrbitState.GetVector();
  ad::Vector6real coeVec = CartToCoe(cartVec, mu);
  ClassicalOE coe(coeVec);
  return coe;
}
ad::Vector6real CartToCoe(const ad::Vector6real &cartVec, double mu) {
  ad::real a, e, p, i, Omega, w, nu, M;

  ad::Vector3real r = cartVec.head(3);
  ad::Vector3real v = cartVec.tail(3);

  ad::real rnorm = norm(r);
  ad::real vnorm = norm(v);

  ad::Vector3real K;
  K << 0, 0, 1.0;

  ad::Vector3real h;
  h = cross(r, v);
  ad::Vector3real n;
  n = cross(K, h);

  ad::real hnorm = norm(h);
  ad::real nnorm = norm(n);

  ad::Vector3real evec =
      ((pow(vnorm, 2.0) - mu / rnorm) * r - dot(r, v) * v) / mu;
  e = norm(evec);

  ad::real xi = pow(vnorm, 2.0) / 2.0 - mu / rnorm;
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

  ad::real ndote = dot(n, evec);
  w = acos(ndote / (nnorm * e));
  if (evec(2).val() < 0) {
    w = 2 * M_PI - w;
  }

  ad::real edotr = dot(evec, r);
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

  ad::Vector6real coeVec;
  coeVec << a, e, i, Omega, w, M;
  return coeVec;
}

/**
 * @brief Compute true anomaly from the eccentric anomaly
 *
 * @param E          Eccentric anomaly [rad]
 * @param e          Eccentricity
 * @return ad::real  True anomaly [rad]
 */
ad::real EccentricAnomToTrueAnom(const ad::real E, const ad::real e) {
  return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e);
  // ad::real beta = e / (1.0 + sqrt(1.0 - pow(e, 2)));
  // return E + 2.0 * atan2(beta * sin(E), 1.0 + beta * cos(E));
}

/**
 * @brief Compute mean anomaly from the eccentric anomaly
 *
 * @param E          Eccentric anomaly [rad]
 * @param e          Eccentricity
 * @return ad::real  Mean anomaly [rad]
 */
ad::real EccentricAnomToMeanAnom(const ad::real E, const ad::real e) {
  return wrapToPi(E - e * sin(E));
}

/**
 * @brief Compute mean anomaly from the Eccentric anomaly
 * @param M          Mean anomaly [rad]
 * @param e          Eccentricity
 * @return ad::real  Eccentric anomaly [rad]
 */
ad::real MeanAnomToEccentricAnom(const ad::real M, const ad::real e) {
  ad::real MM = wrapToPi(M);

  // Initial estimate of E
  ad::real E = MM;
  ad::real Eest = E - (E - e * sin(E) - MM) / (1.0 - e * cos(E));

  double tol = 1e-9;
  int max_itr = 100;
  int itr = 0;

  while ((abs(Eest - E).val() >= tol) && (itr <= max_itr)) {
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
 * @return ad::real  True anomaly [rad]
 */
ad::real MeanAnomToTrueAnom(const ad::real M, const ad::real e) {
  ad::real E = MeanAnomToEccentricAnom(M, e);
  ad::real nu = EccentricAnomToTrueAnom(E, e);

  return wrapToPi(nu);
}

/**
 * @brief Compute true anomaly from the Eccentric anomaly
 * @param nu         True anomaly [rad]
 * @param e          Eccentricity
 * @return ad::real  Eccentric anomaly [rad]
 */
ad::real TrueAnomToEccentricAnom(const ad::real nu, const ad::real e) {
  ad::real E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
  return E;
}

/**
 * @brief Compute mean anomaly from the True anomaly
 *
 * @param nu         True anomaly [rad]
 * @param e          Eccentricity
 * @return ad::real  Mean anomaly [rad]
 */
ad::real TrueAnomToMeanAnom(const ad::real nu, const ad::real e) {
  ad::real E = TrueAnomToEccentricAnom(nu, e);
  ad::real M = EccentricAnomToMeanAnom(E, e);
  return M;
}

ad::Vector6real RoeToCoe(const ad::Vector6real &coe_c,
                         const ad::Vector6real &roe) {
  ad::real ac = coe_c[0];
  ad::real ec = coe_c[1];
  ad::real ic = coe_c[2];
  ad::real Oc = coe_c[3];
  ad::real uc = wrapToPi(coe_c[4] + coe_c[5]);
  ad::real exc = ec * cos(coe_c[4]);
  ad::real eyc = ec * sin(coe_c[4]);

  ad::real dO = roe[5] / (ac * sin(ic));
  ad::real du = (roe[1] / ac) - (dO * cos(ic));

  ad::real ad = ac + roe[0];
  ad::real exd = exc + (roe[2] / ac);
  ad::real eyd = eyc + (roe[3] / ac);
  ad::real id = ic + (roe[4] / ac);
  ad::real Od = Oc + dO;
  ad::real ud = uc + du;

  ad::real wd = atan2(eyd, exd);
  ad::real ed = sqrt(exd * exd + eyd * eyd);
  ad::real Md = wrapToPi(ud - wd);

  ad::Vector6real coe_d;
  coe_d << ad, ed, id, Od, wd, Md;

  return coe_d;
}

ClassicalOE RoeToCoe(const ClassicalOE coe_c, const QuasiNonsingularROE roe) {
  ad::Vector6real coe_d = RoeToCoe(coe_c.GetVector(), roe.GetVector());
  return ClassicalOE(coe_d);
}

ad::Vector6real InertialToRtn(const ad::Vector6real &rtnOrigin,
                              const ad::Vector6real &inertialVec) {
  ad::Vector3real rInertial = inertialVec.head(3);
  ad::Vector3real vInertial = inertialVec.tail(3);
  ad::Vector3real rOrigin = rtnOrigin.head(3);
  ad::Vector3real vOrigin = rtnOrigin.tail(3);

  ad::Vector3real uR, uT, uN;  // RTN basis vectors
  uR = rOrigin.normalized();
  uN = (rOrigin.cross(vOrigin)).normalized();
  uT = uN.cross(uR);

  ad::Matrix3real R_Inertial_Rtn;  // Rotation matrix from inertial to RTN
  R_Inertial_Rtn << uR.transpose(), uT.transpose(), uN.transpose();

  ad::Vector3real w, rRtn, vRtn;
  w = rOrigin.cross(vOrigin) / rOrigin.squaredNorm();
  rRtn = R_Inertial_Rtn * (rInertial - rOrigin);
  vRtn = R_Inertial_Rtn * (vInertial - vOrigin - w.cross(rInertial - rOrigin));

  ad::Vector6real rtnVec;
  rtnVec << rRtn, vRtn;
  return rtnVec;
}

CartesianOrbitState InertialToRtn(
    const CartesianOrbitState &rtnOrigin,
    const CartesianOrbitState &inertialOrbitState) {
  ad::Vector6real rtnOriginVec = rtnOrigin.GetVector();
  ad::Vector6real inertialVec = inertialOrbitState.GetVector();
  ad::Vector6real rtnVec = InertialToRtn(rtnOriginVec, inertialVec);
  CartesianOrbitState rtnOrbitState(rtnVec);
  return rtnOrbitState;
}

ad::Vector6real CoeToRtn(const ad::Vector6real &coe_c,
                         const ad::Vector6real &coe_d, double mu) {
  ad::Vector6real rtnOrigin = CoeToCart(coe_c, mu);
  ad::Vector6real inertialVec = CoeToCart(coe_d, mu);
  ad::Vector6real rtnVec = InertialToRtn(rtnOrigin, inertialVec);
  return rtnVec;
}

CartesianOrbitState CoeToRtn(const ClassicalOE &coe_c, const ClassicalOE &coe_d,
                             double mu) {
  ad::Vector6real rtnVec = CoeToRtn(coe_c.GetVector(), coe_d.GetVector(), mu);
  return CartesianOrbitState(rtnVec);
}

QuasiNonsingularOE CoeToQnsoe(const ClassicalOE &coe) {
  ad::Vector6real qnsoe = CoeToQnsoe(coe.GetVector());
  return QuasiNonsingularOE(qnsoe);
}

ad::Vector6real CoeToQnsoe(const ad::Vector6real &coeVec) {
  ad::real a = coeVec(0);
  ad::real e = coeVec(1);
  ad::real i = coeVec(2);
  ad::real Omega = coeVec(3);
  ad::real w = coeVec(4);
  ad::real M = coeVec(5);

  ad::real u = w + M;
  ad::real ex = e * cos(w);
  ad::real ey = e * sin(w);

  ad::Vector6real qnsoeVec;
  qnsoeVec << a, u, ex, ey, i, Omega;
  return qnsoeVec;
}
ClassicalOE QnsoeToCoe(const QuasiNonsingularOE &qnsoe) {
  ad::Vector6real coeVec = QnsoeToCoe(qnsoe.GetVector());
  ClassicalOE coe(coeVec);
  return coe;
}
ad::Vector6real QnsoeToCoe(const ad::Vector6real &qnsoeVec) {
  ad::real a = qnsoeVec(0);
  ad::real u = qnsoeVec(1);
  ad::real ex = qnsoeVec(2);
  ad::real ey = qnsoeVec(3);
  ad::real i = qnsoeVec(4);
  ad::real Omega = qnsoeVec(5);

  ad::real e = sqrt(ex * ex + ey * ey);
  ad::real w = atan2(ey, ex);
  ad::real M = u - w;

  ad::Vector6real coeVec;
  coeVec << a, e, i, Omega, w, M;
  return coeVec;
}

QuasiNonsingularROE QnsoeToQnsroe(const QuasiNonsingularOE &qnsoe_c,
                                  const QuasiNonsingularOE &qnsoe_d) {
  ad::real a_c = qnsoe_c.a();
  ad::real u_c = qnsoe_c.u();
  ad::real ex_c = qnsoe_c.ex();
  ad::real ey_c = qnsoe_c.ey();
  ad::real i_c = qnsoe_c.i();
  ad::real Omega_c = qnsoe_c.Omega();

  ad::real a_d = qnsoe_d.a();
  ad::real u_d = qnsoe_d.u();
  ad::real ex_d = qnsoe_d.ex();
  ad::real ey_d = qnsoe_d.ey();
  ad::real i_d = qnsoe_d.i();
  ad::real Omega_d = qnsoe_d.Omega();

  ad::real da = (a_d - a_c) / a_c;
  ad::real dl = (u_d - u_c) + (Omega_d - Omega_c) * cos(i_c);
  ad::real dex = ex_d - ex_c;
  ad::real dey = ey_d - ey_c;
  ad::real dix = i_d - i_c;
  ad::real diy = (Omega_d - Omega_c) * sin(i_c);

  ad::Vector6real qnsroe;
  qnsroe << da, dl, dex, dey, dix, diy;
  return QuasiNonsingularROE(qnsroe);
}

ad::Vector6real QnsoeToQnsroe(const ad::Vector6real &qnsoe_c,
                              const ad::Vector6real &qnsoe_d) {
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

ad::Vector6real EquioeToCoe(const ad::Vector6real &equioe) {
  ad::real a = equioe(0);
  ad::real Psi = equioe(1);
  ad::real tq1 = equioe(2);
  ad::real tq2 = equioe(3);
  ad::real p1 = equioe(4);
  ad::real p2 = equioe(5);

  ad::real Omega = atan2(p2, p1);
  ad::real i = 2 * atan2(p1, cos(Omega));

  ad::real wtilde = atan2(tq2, tq1);
  ad::real e = sqrt(tq1 * tq1 + tq2 * tq2);

  ad::real w = wtilde - Omega;
  ad::real f = Psi - wtilde;

  ad::real E = atan2(sin(f) * sqrt(1 - e * e), cos(f) + e);
  ad::real M = E - e * sin(E);

  w = std::fmod(w.val(), 2 * M_PI);
  M = std::fmod(M.val(), 2 * M_PI);
  if (std::abs(M.val() - 2 * M_PI) < std::numeric_limits<double>::epsilon()) {
    M = 0;
  }

  ad::Vector6real coe;
  coe << a, e, i, Omega, w, M;

  return coe;
}

ad::Vector6real MeanToOsculating(const ad::Vector6real &meanCoe, double J2) {
  ad::Vector6real oscCoe;

  if (J2 > 0) {
    ad::Vector6real meanEquioe = CoeToEquioe(meanCoe);
    ad::Vector6real oscEquioe;  // = MeanOscClosedEqui(meanEquioe, J2);
    oscCoe = EquioeToCoe(oscEquioe);
  } else {
    oscCoe = meanCoe;
  }

  return oscCoe;
}

ClassicalOE MeanToOsculating(const ClassicalOE &meanCoe, double J2) {
  ad::Vector6real meanCoeVec = meanCoe.GetVector();
  ad::Vector6real oscCoeVec = MeanToOsculating(meanCoeVec, J2);
  ClassicalOE oscCoe(oscCoeVec);
  return oscCoe;
}

ad::Vector6real osc2mean_NRiterator(const ad::Vector6real &osc_equi_elem,
                                    double tol) {
  ad::Vector6real mean_equi_elem = osc_equi_elem;
  double R = 1.0;
  int niter = 0;

  while (std::abs(R) > tol) {
    niter++;
    ad::Vector6real osc_loop;
    // std::tie(std::ignore, osc_loop, std::ignore) =
    // transformationmatrix_osc2mean_equinoctial(mean_equi_elem);
    ad::Vector6real delta = osc_equi_elem - osc_loop;
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

ad::Vector6real OsculatingToMean(const ad::Vector6real &oscCoe, double J2) {
  ad::Vector6real meanCoe;
  double tol = 1e-8;

  if (J2 > 0) {
    ad::Vector6real oscEquioe = CoeToEquioe(oscCoe);
    ad::Vector6real meanEquioe = osc2mean_NRiterator(oscEquioe, tol);
    meanCoe = EquioeToCoe(meanEquioe);
  } else {
    meanCoe = oscCoe;
  }

  return meanCoe;
}

ClassicalOE OsculatingToMean(const ClassicalOE &oscCoe, double J2) {
  ad::Vector6real oscCoeVec = oscCoe.GetVector();
  ad::Vector6real meanCoeVec = OsculatingToMean(oscCoeVec, J2);
  ClassicalOE meanCoe(meanCoeVec);
  return meanCoe;
}

ad::Vector6real CoeToEquioe(const ad::Vector6real &coe) {
  ad::real a = coe(0);
  ad::real e = coe(1);
  ad::real i = coe(2);
  ad::real Omega = coe(3);
  ad::real w = coe(4);
  ad::real M = coe(5);

  ad::real f = MeanAnomToTrueAnom(M, e);
  ad::real w_tilde = Omega + w;
  ad::real Psi = w_tilde + f;
  ad::real tq1 = e * cos(w_tilde);
  ad::real tq2 = e * sin(w_tilde);
  ad::real p1 = tan(i / 2) * cos(Omega);
  ad::real p2 = tan(i / 2) * sin(Omega);

  Psi = fmod(Psi.val(), 2 * M_PI);
  if (Psi > M_PI) {
    Psi = Psi - 2 * M_PI;
  }

  ad::Vector6real equioe;
  equioe << a, Psi, tq1, tq2, p1, p2;

  return equioe;
}

ad::Vector6real CartToQnsoe(const ad::Vector6real &cart, double mu) {
  ad::Vector6real coe = CartToCoe(cart, mu);
  ad::Vector6real qnsoe = CoeToQnsoe(coe);
  return qnsoe;
}
QuasiNonsingularOE CartToQnsoe(const CartesianOrbitState &cart, double mu) {
  ad::Vector6real cartVec = cart.GetVector();
  ad::Vector6real qnsoeVec = CartToQnsoe(cartVec, mu);
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

ad::Vector6real DelaunayToCoe(const ad::Vector6real &delaunay, double mu,
                              double n, double t) {
  ad::real l = delaunay(0);
  ad::real g = delaunay(1);
  ad::real h = delaunay(2);
  ad::real L = delaunay(3);
  ad::real G = delaunay(4);
  ad::real H = delaunay(5);

  ad::real a = L * L / mu;
  ad::real e = sqrt(1 - pow(G / L, 2));
  ad::real i = acos(H / G);
  ad::real O = h + n * t;
  ad::real w = g;
  ad::real M = l;

  ad::Vector6real coe;
  coe << a, e, i, O, w, M;
  return coe;
}

ad::Vector6real CoeToDelaunay(const ad::Vector6real &coe, double mu, double n,
                              double t) {
  ad::real a = coe(0);
  ad::real e = coe(1);
  ad::real i = coe(2);
  ad::real O = coe(3);
  ad::real w = coe(4);
  ad::real M = coe(5);

  ad::real l = M;
  ad::real g = w;
  ad::real h = O - n * t;
  ad::real L = sqrt(mu * a);
  ad::real G = L * sqrt(1 - e * e);
  ad::real H = G * cos(i);

  ad::Vector6real delaunay;
  delaunay << l, g, h, L, G, H;
  return delaunay;
}

}  // namespace LPT