#include "lupnt/measurements/surface_measurements.h"

#include "lupnt/core/constants.h"

namespace lupnt {

  double SurfaceLansMeasurement::PredictedRange(const Vec3d& r_rover, double clock_bias_s) const {
    double range = (r_sat - r_rover).norm();
    return range + C * clock_bias_s;
  }

  Vec3d SurfaceLansMeasurement::LosUnit(const Vec3d& r_rover) const {
    Vec3d los = r_sat - r_rover;
    double range = los.norm();
    return (range > 0.0) ? Vec3d(los / range) : Vec3d(Vec3d::Zero());
  }

  MeasData SurfaceLansMeasurement::Compute(const NavErrorContext& nom, MatXd* H) const {
    MeasData md;
    md.timestamp = timestamp;
    md.value = VecXd::Constant(1, PredictedRange(nom.r, nom.clock_bias_s));
    md.covariance = MatXd::Constant(1, 1, NoiseVariance());
    if (H != nullptr) {
      Vec3d u = LosUnit(nom.r);  // rover -> satellite unit vector
      *H = MatXd::Zero(1, nom.error_state_size);
      H->block(0, config.i_dr, 1, 3) = -u.transpose();  // d(range)/d(r)
      (*H)(0, config.i_dcb) = C;                        // d(rho)/d(clock_bias)
    }
    return md;
  }

}  // namespace lupnt
