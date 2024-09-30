#include "lupnt/lupnt.h"

using namespace lupnt;

constexpr double ABS_TOL = 1e-6;
constexpr double REL_TOL = 1e-6;

static Vec6 GetClassicalOE() {
  Real a = 5740;        // [km] Semi-major axis
  Real e = 0.58;        // [-] Eccentricity
  Real i = 54.9 * RAD;  // [rad] Inclination
  Real O = 0;           // [rad] Right Ascension of Ascending Node
  Real w = 86.3 * RAD;  // [rad] Argument of Perigee
  Real M = 0;           // [rad] True Anomaly

  Vec6 coe = {a, e, i, O, w, M};
  return coe;
}

static MatX6 GetClassicalOEMat(int n) {
  Vec6 coe = GetClassicalOE();
  MatX6 coe_mat = coe.transpose().replicate(n, 1);
  coe_mat.col(5) = VecX::LinSpaced(n, 0, TWO_PI);
  return coe_mat;
}
