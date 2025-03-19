#include <cassert>
#include <iostream>
#include <vector>

#include "lupnt/core/constants.h"

namespace lupnt {

  double LinearInterp1d(const VecXd& x, const VecXd& data, double ix);
  double LinearInterp2d(const VecXd& x, const VecXd& y, const MatXd& data, double ix, double iy);

  class LagrangeInterpolator {
  public:
    LagrangeInterpolator(const VecXd& x, double xi, int order);
    double Interpolate(const VecXd& data);

  private:
    VecXd x_;
    double xi_;
    int order_;
    VecXd weights_;
    int i0_;

    void ComputeFirstIndex();
    void ComputeWeights();
  };

}  // namespace lupnt
