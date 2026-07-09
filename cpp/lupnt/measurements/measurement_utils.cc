#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

#include "lupnt/core/constants.h"

namespace lupnt {

  Vec1 Range(const VecX& r1, const VecX& r2) { return Vec1((r2 - r1).norm()); }

  Vec1 RangeRate(const VecX& r1, const VecX& r2, const VecX& v1, const VecX& v2) {
    return Vec1((v2 - v1).dot(r2 - r1) / (r2 - r1).norm());
  }

  Vec1 RangeRate(const VecX& rv1, const VecX& rv2) {
    const int dim = rv1.size() / 2;
    VecX r1 = rv1.head(dim), r2 = rv2.head(dim);
    VecX v1 = rv1.tail(dim), v2 = rv2.tail(dim);
    return RangeRate(r1, r2, v1, v2);
  }

  Vec2 RangeAndRangeRate(const VecX& r1, const VecX& r2, const VecX& v1, const VecX& v2) {
    Real r_norm = (r2 - r1).norm();
    return Vec2(r_norm, (v2 - v1).dot(r2 - r1) / r_norm);
  }

  Vec2 RangeAndRangeRate(const VecX& rv1, const VecX& rv2) {
    const int dim = rv1.size() / 2;
    VecX r1 = rv1.head(dim), r2 = rv2.head(dim);
    VecX v1 = rv1.tail(dim), v2 = rv2.tail(dim);
    return RangeAndRangeRate(r1, r2, v1, v2);
  }

  Vec1 Pseudorange(const VecX& r1, const VecX& r2, Real b1, Real b2) {
    return Range(r1, r2).array() + C * (b1 - b2);
  }

  Vec1 PseudorangeRate(const VecX& r1, const VecX& r2, const VecX& v1, const VecX& v2, Real d1,
                       Real d2) {
    return RangeRate(r1, r2, v1, v2).array() + C * (d1 - d2);
  }

  Vec1 PseudorangeRate(const VecX& rv1, const VecX& rv2, Real d1, Real d2) {
    const int dim = rv1.size() / 2;
    VecX r1 = rv1.head(dim), r2 = rv2.head(dim);
    VecX v1 = rv1.tail(dim), v2 = rv2.tail(dim);
    return PseudorangeRate(r1, r2, v1, v2, d1, d2);
  }

}  // namespace lupnt
