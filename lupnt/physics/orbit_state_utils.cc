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

#include <functional>
#include <map>
#include <queue>
#include <tuple>

#include "lupnt/numerics/math_utils.h"

#define ABSOLUTE_CONVERSION(from, to, func)                                \
  {                                                                        \
    {OrbitStateRepres::from, OrbitStateRepres::to},                        \
        [](const Vector6 &x, double mu) -> Vector6 { return func(x, mu); } \
  }

#define RELATIVE_CONVERSION(from, to, func)                 \
  {                                                         \
    {OrbitStateRepres::from, OrbitStateRepres::to},         \
        [](const Vector6 &x, const Vector6 &y) -> Vector6 { \
          return func(x, y);                                \
        }                                                   \
  }

namespace lupnt {

namespace {

std::vector<OrbitStateRepres> FindShortestPath(
    OrbitStateRepres start, OrbitStateRepres end,
    const std::map<std::tuple<OrbitStateRepres, OrbitStateRepres>,
                   std::function<Vector6(const Vector6 &, double)>>
        &absolute_conversions) {
  std::queue<OrbitStateRepres> queue;
  std::map<OrbitStateRepres, OrbitStateRepres> predecessors;
  std::map<OrbitStateRepres, bool> visited;

  queue.push(start);
  visited[start] = true;
  predecessors[start] = start;  // Start node is its own predecessor

  while (!queue.empty()) {
    OrbitStateRepres current = queue.front();
    queue.pop();

    if (current == end) {
      // Path found, reconstruct it
      std::vector<OrbitStateRepres> path;
      for (OrbitStateRepres at = end; at != start; at = predecessors[at]) {
        path.push_back(at);
      }
      path.push_back(start);
      std::reverse(path.begin(), path.end());
      return path;
    }

    // Explore neighbors
    for (const auto &entry : absolute_conversions) {
      OrbitStateRepres neighbor = std::get<1>(entry.first);
      if (std::get<0>(entry.first) == current && !visited[neighbor]) {
        queue.push(neighbor);
        visited[neighbor] = true;
        predecessors[neighbor] = current;
      }
    }
  }

  throw std::runtime_error("Path not found from start to end representation.");
}

}  // namespace

std::map<std::tuple<OrbitStateRepres, OrbitStateRepres>,
         std::function<Vector6(const Vector6 &, double)>>
    absolute_conversions = {
        ABSOLUTE_CONVERSION(CARTESIAN, CLASSICAL_OE, CartesianToClassical),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, CARTESIAN, ClassicalToCartesian),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, QUASI_NONSINGULAR_OE,
                            ClassicalToQuasiNonsingular),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, EQUINOCTIAL_OE,
                            ClassicalToEquinoctial),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, DELAUNAY_OE, ClassicalToDelaunay),
        ABSOLUTE_CONVERSION(QUASI_NONSINGULAR_OE, CLASSICAL_OE,
                            QuasiNonsingularToClassical),
        ABSOLUTE_CONVERSION(EQUINOCTIAL_OE, CLASSICAL_OE,
                            EquinoctialToClassical),
        ABSOLUTE_CONVERSION(DELAUNAY_OE, CLASSICAL_OE, DelaunayToClassical),
};

std::map<std::tuple<OrbitStateRepres, OrbitStateRepres>,
         std::function<Vector6(const Vector6 &, const Vector6 &)>>
    relative_conversions = {
        RELATIVE_CONVERSION(CARTESIAN, RTN, InertialToRtn),
        RELATIVE_CONVERSION(QUASINONSINGULAR_ROE, CLASSICAL_OE,
                            RelativeQuasiNonsingularToClassical),
};

Vector6 ConvertOrbitState(const Vector6 &state_in, OrbitStateRepres repres_in,
                          OrbitStateRepres repres_out, double mu) {
  if (repres_in == repres_out) {
    return state_in;
  }

  std::vector<OrbitStateRepres> path =
      FindShortestPath(repres_in, repres_out, absolute_conversions);

  Vector6 state = state_in;
  for (size_t i = 0; i < path.size() - 1; i++) {
    state = absolute_conversions[{path[i], path[i + 1]}](state, mu);
  }

  return state;
}

Vector6 ConvertOrbitState(const Vector6 &state_in_c, const Vector6 &state_in_d,
                          OrbitStateRepres repres_in_c,
                          OrbitStateRepres repres_in_d,
                          OrbitStateRepres repres_out, double mu) {
  // Check case
  // - (absolute_c, absolute_d) to relative_d
  // - (absolute_c, relative_d) to absolute_d
  bool to_relative =
      repres_in_c < OrbitStateRepres::ABSOLUTE_RELATIVE_SEPARATOR &&
      repres_in_d < OrbitStateRepres::ABSOLUTE_RELATIVE_SEPARATOR;

  for (const auto &entry : relative_conversions) {
    auto repres_from = std::get<0>(entry.first);
    auto repres_to = std::get<1>(entry.first);
    auto func = entry.second;

    if (to_relative && repres_to == repres_out) {
      // state_in_c is absolute
      // state_in_d is absolute
      auto state_abs_c =
          ConvertOrbitState(state_in_c, repres_in_c, repres_from, mu);
      auto state_abs_d =
          ConvertOrbitState(state_in_d, repres_in_d, repres_from, mu);
      return func(state_abs_c, state_abs_d);
    } else if (!to_relative && repres_from == repres_in_d) {
      // state_in_c is absolute
      // state_in_d is relative
      auto state_abs_c =
          ConvertOrbitState(state_in_c, repres_in_c, repres_to, mu);
      auto state_abs_d = func(state_abs_c, state_in_d);
      return ConvertOrbitState(state_abs_d, repres_to, repres_out, mu);
    }
  }
  throw std::runtime_error(
      "Relative conversion not found for the given input.");
}

std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
    const std::shared_ptr<OrbitState> &state_in, OrbitStateRepres repres_out,
    double mu) {
  Vector6 state_out = ConvertOrbitState(
      state_in->GetVector(), state_in->GetOrbitStateRepres(), repres_out, mu);
  return std::make_shared<OrbitState>(state_out, state_in->GetCoordSystem(),
                                      repres_out, state_in->GetNames(),
                                      state_in->GetUnits());
}

// From CartesianOrbitState
// - To ClassicalOE
ClassicalOE CartesianToClassical(const CartesianOrbitState &rv, double mu) {
  return ClassicalOE(CartesianToClassical(rv.GetVector(), mu),
                     rv.GetCoordSystem());
}

Vector6 CartesianToClassical(const Vector6 &rv, double mu) {
  Vector3 r = rv.head(3);
  Vector3 v = rv.tail(3);

  Vector3 k{0, 0, 1};

  Vector3 h = r.cross(v);
  Vector3 n = k.cross(h);

  Vector3 evec = v.cross(h) / mu - r / r.norm();

  real e = evec.norm();
  real i = safe_acos(h(2) / h.norm());
  real a = 1.0 / (2.0 / r.norm() - pow(v.norm(), 2.0) / mu);

  real nu = (r.dot(v) >= 0) ? angleBetweenVectors(evec, r)
                            : -angleBetweenVectors(evec, r);

  real Omega =
      (n(1) >= 0) ? safe_acos(n(0) / n.norm()) : -safe_acos(n(0) / n.norm());

  real w = (evec(2) >= 0) ? angleBetweenVectors(n, evec)
                          : -angleBetweenVectors(n, evec);

  real M = TrueToMeanAnomaly(nu, e);

  return Vector6{a, e, i, Omega, w, M};
}

// - To CartesianOrbitState (relative)
CartesianOrbitState InertialToRtn(const CartesianOrbitState &rv_c,
                                  const CartesianOrbitState &rv_d) {
  return CartesianOrbitState(InertialToRtn(rv_c.GetVector(), rv_d.GetVector()),
                             rv_c.GetCoordSystem());
}

Vector6 InertialToRtn(const Vector6 &rv_c, const Vector6 &rv_d) {
  Vector3 r_d = rv_d.head(3);
  Vector3 v_d = rv_d.tail(3);
  Vector3 r_c = rv_c.head(3);
  Vector3 v_c = rv_c.tail(3);

  // RTN basis vectors
  Vector3 uR = r_c.normalized();
  Vector3 uN = (r_c.cross(v_c)).normalized();
  Vector3 uT = uN.cross(uR);

  Matrix3 Rot_inert_rtn;  // Rotation matrix from inertial to RTN
  Rot_inert_rtn << uR.transpose(), uT.transpose(), uN.transpose();

  Vector3 w = r_c.cross(v_c) / r_c.norm();
  Vector3 r_rtn_d = Rot_inert_rtn * (r_d - r_c);
  Vector3 v_rtn_d = Rot_inert_rtn * (v_d - v_c - w.cross(r_d - r_c));

  Vector6 rv_rtn_d;
  rv_rtn_d << r_rtn_d, v_rtn_d;
  return rv_rtn_d;
}

CartesianOrbitState RtnToInertial(const CartesianOrbitState &rv_c,
                                  const CartesianOrbitState &rv_rtn_d) {
  return CartesianOrbitState(
      RtnToInertial(rv_c.GetVector(), rv_rtn_d.GetVector()),
      rv_c.GetCoordSystem());
}

Vector6 RtnToInertial(const Vector6 &rv_c, const Vector6 &rv_rtn_d) {
  Vector3 r_c = rv_c.head(3);
  Vector3 v_c = rv_c.tail(3);

  // RTN basis vectors
  Vector3 uR = r_c.normalized();
  Vector3 uN = (r_c.cross(v_c)).normalized();
  Vector3 uT = uN.cross(uR);

  Matrix3 Rot_rtn_inert;  // Rotation matrix from RTN to inertial
  Rot_rtn_inert << uR, uT, uN;

  Vector3 w = r_c.cross(v_c) / r_c.norm();
  Vector3 r_d = r_c + Rot_rtn_inert * rv_rtn_d.head(3);
  Vector3 v_d = v_c + Rot_rtn_inert * rv_rtn_d.tail(3) + w.cross(r_d - r_c);

  Vector6 rv_d;
  rv_d << r_d, v_d;
  return rv_d;
}

// From ClassicalOE
// - To CartesianOrbitState
CartesianOrbitState ClassicalToCartesian(const ClassicalOE &coe, double mu) {
  return CartesianOrbitState(ClassicalToCartesian(coe.GetVector(), mu),
                             coe.GetCoordSystem());
}

Vector6 ClassicalToCartesian(const Vector6 &coe, double mu) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  real p = a * (1.0 - pow(e, 2.0));
  real nu = MeanToTrueAnomaly(M, e);
  real pev = p / (1.0 + e * cos(nu));
  real mu_p = sqrt(mu / p);

  Vector3 r_PQW = {pev * cos(nu), pev * sin(nu), 0.0};
  Vector3 v_PQW = {-mu_p * sin(nu), mu_p * (e + cos(nu)), 0.0};

  // rot = Rot3(-Omega) * Rot1(-i) * Rot3(-w)
  Matrix3 rot{
      {cos(Omega) * cos(w) - sin(Omega) * sin(w) * cos(i),
       -cos(Omega) * sin(w) - sin(Omega) * cos(w) * cos(i),
       sin(Omega) * sin(i)},
      {sin(Omega) * cos(w) + cos(Omega) * sin(w) * cos(i),
       -sin(Omega) * sin(w) + cos(Omega) * cos(w) * cos(i),
       -cos(Omega) * sin(i)},
      {sin(w) * sin(i), cos(w) * sin(i), cos(i)},
  };

  Vector3 r = rot * r_PQW;
  Vector3 v = rot * v_PQW;

  Vector6 rv;
  rv << r, v;
  return rv;
}

// - To QuasiNonsingularOE
QuasiNonsingularOE ClassicalToQuasiNonsingular(const ClassicalOE &coe,
                                               double mu) {
  return QuasiNonsingularOE(ClassicalToQuasiNonsingular(coe.GetVector()),
                            coe.GetCoordSystem());
}

Vector6 ClassicalToQuasiNonsingular(const Vector6 &coe, double mu) {
  auto [a, e, i, Omega, w, M] = unpack(coe);

  real u = w + M;
  real ex = e * cos(w);
  real ey = e * sin(w);

  return Vector6(a, u, ex, ey, i, Omega);
}

// - To EquinoctialOE
EquinoctialOE ClassicalToEquinoctial(const ClassicalOE &coe, double mu) {
  return EquinoctialOE(ClassicalToEquinoctial(coe.GetVector()),
                       coe.GetCoordSystem());
}

Vector6 ClassicalToEquinoctial(const Vector6 &coe, double mu) {
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

  return Vector6(a, Psi, tq1, tq2, p1, p2);
}

// - To DelaunayOE
DelaunayOE ClassicalToDelaunay(const ClassicalOE &coe, double mu) {
  return DelaunayOE(ClassicalToDelaunay(coe.GetVector(), mu),
                    coe.GetCoordSystem());
}

Vector6 ClassicalToDelaunay(const Vector6 &coe, double mu) {
  auto [a, e, i, O, w, M] = unpack(coe);

  real n = sqrt(mu / pow(a, 3));
  real t = M / n;

  real l = M;
  real g = w;
  real h = O - n * t;
  real L = sqrt(mu * a);
  real G = L * sqrt(1 - e * e);
  real H = G * cos(i);

  return Vector6(l, g, h, L, G, H);
}

// From QuasiNonsingularOE
// - To ClassicalOE
ClassicalOE QuasiNonsingularToClassical(const QuasiNonsingularOE &qnsoe,
                                        double mu) {
  return ClassicalOE(QuasiNonsingularToClassical(qnsoe.GetVector()),
                     qnsoe.GetCoordSystem());
}

Vector6 QuasiNonsingularToClassical(const Vector6 &qnsoeVec, double mu) {
  auto [a, u, ex, ey, i, Omega] = unpack(qnsoeVec);

  real e = sqrt(ex * ex + ey * ey);
  real w = atan2(ey, ex);
  real M = u - w;

  Vector6 coe{a, e, i, Omega, w, M};
  return coe;
}

// From EquinoctialOE
// - To ClassicalOE
ClassicalOE EquinoctialToClassical(const EquinoctialOE &eqoe, double mu) {
  return ClassicalOE(EquinoctialToClassical(eqoe.GetVector()),
                     eqoe.GetCoordSystem());
}

Vector6 EquinoctialToClassical(const Vector6 &equioe, double mu) {
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

  return Vector6(a, e, i, Omega, w, M);
}

// From DelaunayOE
// - To ClassicalOE
ClassicalOE DelaunayToClassical(const DelaunayOE &deloe, double mu) {
  return ClassicalOE(DelaunayToClassical(deloe.GetVector(), mu),
                     deloe.GetCoordSystem());
}

Vector6 DelaunayToClassical(const Vector6 &delaunay, double mu) {
  auto [l, g, h, L, G, H] = unpack(delaunay);

  real a = L * L / mu;
  real M = l;

  real n = sqrt(mu / pow(a, 3));
  real t = M / n;

  real e = sqrt(1 - pow(G / L, 2));
  real i = safe_acos(H / G);
  real O = h + n * t;
  real w = g;

  return Vector6(a, e, i, O, w, M);
}

// From ClassicalOE QuasiNonsingularOE
// - To ClassicalOE
ClassicalOE RelativeQuasiNonsingularToClassical(
    const ClassicalOE &coe_c,
    const QuasiNonsingularROE &RelativeQuasiNonsingular) {
  return ClassicalOE(
      RelativeQuasiNonsingularToClassical(coe_c.GetVector(),
                                          RelativeQuasiNonsingular.GetVector()),
      coe_c.GetCoordSystem());
}

Vector6 RelativeQuasiNonsingularToClassical(
    const Vector6 &coe_c, const Vector6 &RelativeQuasiNonsingular) {
  auto [ac, ec, ic, Oc, wc, Mc] = unpack(coe_c);
  auto [ada, adl, adex, adey, adix, adiy] = unpack(RelativeQuasiNonsingular);

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

  return Vector6(ad, ed, id, Od, wd, Md);
}

// Mean and Osculating
Vector6 MeanToOsculating(const Vector6 &coe_m, double J2) {
  Vector6 coe_o;

  if (J2 > 0) {
    Vector6 meanEquioe = ClassicalToEquinoctial(coe_m);
    Vector6 oscEquioe;  // = MeanOscClosedEqui(meanEquioe, J2);
    coe_o = EquinoctialToClassical(oscEquioe);
  } else {
    coe_o = coe_m;
  }

  return coe_o;
}

ClassicalOE MeanToOsculating(const ClassicalOE &coe_m, double J2) {
  return ClassicalOE(MeanToOsculating(coe_m.GetVector(), J2),
                     coe_m.GetCoordSystem());
}

Vector6 osc2mean_NRiterator(const Vector6 &osc_equi_elem, double tol) {
  Vector6 mean_equi_elem = osc_equi_elem;
  double R = 1.0;
  int niter = 0;

  while (std::abs(R) > tol) {
    niter++;
    Vector6 osc_loop;
    // std::tie(std::ignore, osc_loop, std::ignore) =
    // transformationmatrix_osc2mean_equinoctial(mean_equi_elem);
    Vector6 delta = osc_equi_elem - osc_loop;
    // R = delta.norm_inf(); // Assuming the library provides an infinity
    // norm function
    mean_equi_elem = mean_equi_elem + delta;

    if (niter > 100) {
      std::cout << "Osc2Mean iterations > 100" << std::endl;
      break;
    }
  }

  return mean_equi_elem;
}

Vector6 OsculatingToMean(const Vector6 &coe_o, double J2) {
  Vector6 coe_m;
  double tol = 1e-8;

  if (J2 > 0) {
    Vector6 eqoe_o = ClassicalToEquinoctial(coe_o);
    Vector6 eqoe_m = osc2mean_NRiterator(eqoe_o, tol);
    coe_m = EquinoctialToClassical(eqoe_m);
  } else {
    coe_m = coe_o;
  }

  return coe_m;
}

ClassicalOE OsculatingToMean(const ClassicalOE &coe_o, double J2) {
  return ClassicalOE(OsculatingToMean(coe_o.GetVector(), J2),
                     coe_o.GetCoordSystem());
}

real EccentricToTrueAnomaly(real E, real e) {
  return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e);
}

// Anomaly
real EccentricToMeanAnomaly(real E, real e) { return wrapToPi(E - e * sin(E)); }

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

real MeanToTrueAnomaly(real M, real e) {
  real E = MeanToEccentricAnomaly(M, e);
  real nu = EccentricToTrueAnomaly(E, e);
  return wrapToPi(nu);
}

real TrueToEccentricAnomaly(real nu, real e) {
  real E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
  return E;
}

real TrueToMeanAnomaly(real nu, real e) {
  real E = TrueToEccentricAnomaly(nu, e);
  real M = EccentricToMeanAnomaly(E, e);
  return M;
}

}  // namespace lupnt