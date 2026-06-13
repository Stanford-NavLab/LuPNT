#include "lupnt/dynamics/imu_dynamics.h"

#include "lupnt/core/constants.h"
#include "lupnt/devices/imu.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/states/state.h"

namespace lupnt {

  ImuDynamics::ImuDynamics() { imu_model_ = ImuModel::UNDEFINED; }

  ImuDynamics::ImuDynamics(Config& config) {
    imu_model_ = ImuModel::UNDEFINED;
    if (config["model"]) {
      imu_model_ = enum_cast<ImuModel>(config["model"].as<std::string>()).value();
    }
  }

  void ImuDynamics::SetModel(ImuModel imu_model) { imu_model_ = imu_model; }
  ImuModel ImuDynamics::GetModel() const { return imu_model_; }

  State ImuDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    LUPNT_CHECK(x0.GetType() == ImuState::TYPE || x0.GetType() == UNDEFINED,
                "State must be of type ImuState", "ImuDynamics");
    LUPNT_CHECK(x0.size() == 6, "State must be 6-dimensional", "ImuDynamics");
    LUPNT_CHECK(u == nullptr, "ImuDynamics does not support control inputs", "ImuDynamics");

    ImuParameters imu_params = Imu::GetParameters(imu_model_);
    Real sqrt_f = 1.0 / sqrt(tf - t0);
    ImuState imu_f;
    ImuState imu_0 = x0;
    for (int i = 0; i < 3; i++) {
      imu_f.b_w()(i) += SampleNormal(imu_0.b_w()(i), imu_params.sigma_w_bias / sqrt_f, rng_.get());
      imu_f.b_a()(i) += SampleNormal(imu_0.b_a()(i), imu_params.sigma_a_bias / sqrt_f, rng_.get());
    }
    return imu_f;
  }
}  // namespace lupnt
