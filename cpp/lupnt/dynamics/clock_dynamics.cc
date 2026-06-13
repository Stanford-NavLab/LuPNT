#include "lupnt/dynamics/clock_dynamics.h"

#include <tuple>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/devices/clock.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  // ClockDynamics
  // A. C. Long, ‘Goddard Enhanced Onboard Navigation System (GEONS) Mathematical
  // Specifications’, 2024.
  std::tuple<double, double, double> ClockDynamics::GetClockValues(ClockModel clk_model) {
    double q1 = 0.0;  // [s]
    double q2 = 0.0;  // [1/s]
    double q3 = 0.0;  // [1/s^2]

    switch (clk_model) {
      case ClockModel::OCXO:
        // Preliminary Navigation System Design for the First LCRNS
        // Satellite Providing Lunar PNT Services
        q1 = pow(6.2299445014e-13, 2);
        q2 = pow(2.0129544799e-14, 2);
        q3 = pow(7.0118586804e-28, 2);
        break;
      case ClockModel::USO:
        // MMS USO
        q1 = 1.2e-22;
        q2 = 1.58e-26;
        q3 = 5.9e-75;
        break;
      case ClockModel::CSAC:
        q1 = 5.8e-22;
        q2 = 1.5e-29;
        q3 = 1e-45;
        break;
      case ClockModel::MINI_RAFS:
        q1 = pow(1.0e-11, 2);
        q2 = pow(1.1e-15, 2);
        q3 = pow(7.0118586804e-28, 2);
        break;
      case ClockModel::RAFS:
        // https://www.nasa.gov/wp-content/uploads/2025/08/gps-based-autonomous-navigation-study-for-the-lunar-gateway-1.pdf?emrc=68cc48e722112
        q1 = 3.70e-24;
        q2 = 1.87e-33;
        q3 = 7.56e-59;
        break;
      case ClockModel::DSAC:
        q1 = 4.23e-26;
        q2 = 6.19e-38;
        q3 = 2.24e-61;
        break;
      default: LUPNT_CHECK(false, "Unknown clock model", "ClockDynamics");
    }
    return std::make_tuple(q1, q2, q3);
  }

  ClockDynamics::ClockDynamics() { model_ = ClockModel::UNDEFINED; }

  ClockDynamics::ClockDynamics(Config& config) {
    model_ = ClockModel::UNDEFINED;
    if (config["model"]) {
      model_ = enum_cast<ClockModel>(config["model"].as<std::string>()).value();
    }
    if (config["clock_bias_unit"]) {
      bias_unit_ = enum_cast<ClockBiasUnit>(config["clock_bias_unit"].as<std::string>()).value();
    } else if (config["bias_unit"]) {
      bias_unit_ = enum_cast<ClockBiasUnit>(config["bias_unit"].as<std::string>()).value();
    }
  }

  double ClockDynamics::SecondsToBiasUnitScale(ClockBiasUnit unit) {
    switch (unit) {
      case ClockBiasUnit::SECONDS: return 1.0;
      case ClockBiasUnit::METERS: return C;
      case ClockBiasUnit::KILOMETERS: return C * KM_M;
      default: LUPNT_CHECK(false, "Unknown clock bias unit", "ClockDynamics");
    }
  }

  Real ClockDynamics::SecondsToBiasUnits(Real value_s, ClockBiasUnit unit) {
    return value_s * SecondsToBiasUnitScale(unit);
  }

  Real ClockDynamics::BiasUnitsToSeconds(Real value, ClockBiasUnit unit) {
    return value / SecondsToBiasUnitScale(unit);
  }

  std::vector<std::string> ClockDynamics::GetStateUnits(int state_size, ClockBiasUnit unit) {
    std::vector<std::string> units;
    switch (unit) {
      case ClockBiasUnit::SECONDS: units = {"s", "s/s", "s/s^2"}; break;
      case ClockBiasUnit::METERS: units = {"m", "m/s", "m/s^2"}; break;
      case ClockBiasUnit::KILOMETERS: units = {"km", "km/s", "km/s^2"}; break;
      default: LUPNT_CHECK(false, "Unknown clock bias unit", "ClockDynamics");
    }

    LUPNT_CHECK(state_size == 2 || state_size == 3, "Clock state size must be 2 or 3",
                "ClockDynamics");
    units.resize(state_size);
    return units;
  }

  Real ClockDynamics::RelativisticRateCorrection(const Vec3& r_centered, const Vec3& v_centered,
                                                 Real GM, Real c, Real reference_rate_offset) {
    Real r = r_centered.norm();
    LUPNT_CHECK(r > 0.0, "Relativistic clock correction requires a nonzero radius",
                "ClockDynamics");

    Real v2 = v_centered.squaredNorm();
    return reference_rate_offset - (GM / r + 0.5 * v2) / (c * c);
  }

  Mat2 ClockDynamics::TwoStateNoise(ClockModel clk_model, Real dt) {
    Real q1, q2, q3;
    std::tie(q1, q2, q3) = GetClockValues(clk_model);

    Real dt2 = dt * dt;
    Real dt3 = dt2 * dt;

    Mat2 Q;
    Q(0, 0) = q1 * dt + q2 * dt3 / 3.0;
    Q(0, 1) = q2 * dt2 / 2.0;
    Q(1, 0) = Q(0, 1);
    Q(1, 1) = q2 * dt;
    return Q;
  };

  Mat2 ClockDynamics::TwoStateNoise(ClockModel clk_model, Real dt, ClockBiasUnit unit) {
    Real scale = SecondsToBiasUnitScale(unit);
    return TwoStateNoise(clk_model, dt) * scale * scale;
  }

  Mat3 ClockDynamics::ThreeStateNoise(ClockModel clk_model, Real dt) {
    Real q1, q2, q3;
    std::tie(q1, q2, q3) = GetClockValues(clk_model);

    Real dt2 = dt * dt;
    Real dt3 = dt2 * dt;
    Real dt4 = dt3 * dt;
    Real dt5 = dt4 * dt;

    Mat3 Q;
    Q(0, 0) = q1 * dt + q2 * dt3 / 3 + q3 * dt5 / 20.0;
    Q(0, 1) = q2 * dt2 / 2.0 + q3 * dt4 / 8.0;
    Q(0, 2) = q3 * dt3 / 6.0;
    Q(1, 0) = Q(0, 1);
    Q(1, 1) = q2 * dt + q3 * dt3 / 3.0;
    Q(1, 2) = q3 * dt2 / 2.0;
    Q(2, 0) = Q(0, 2);
    Q(2, 1) = Q(1, 2);
    Q(2, 2) = q3 * dt;
    return Q;
  };

  Mat3 ClockDynamics::ThreeStateNoise(ClockModel clk_model, Real dt, ClockBiasUnit unit) {
    Real scale = SecondsToBiasUnitScale(unit);
    return ThreeStateNoise(clk_model, dt) * scale * scale;
  }

  VecX ClockDynamics::GetProcessNoise(Real dt, int state_size) {
    MatX Q;
    VecX noise(state_size);
    if (state_size == 2) {
      Q = ClockDynamics::TwoStateNoise(model_, dt, bias_unit_);
      noise = SampleMvNormal(Vec2d::Zero(), Q, 1, rng_.get()).transpose();
    } else if (state_size == 3) {
      Q = ClockDynamics::ThreeStateNoise(model_, dt, bias_unit_);
      noise = SampleMvNormal(Vec3d::Zero(), Q, 1, rng_.get()).transpose();
    } else {
      LUPNT_CHECK(false, fmt::format("Invalid clock state size {}", state_size), "ClockDynamics");
    }
    return noise;
  }

  Mat2 ClockDynamics::TwoStatePhi(Real dt) {
    Mat2 Phi_clk{{1, dt}, {0, 1}};
    return Phi_clk;
  }

  Mat3 ClockDynamics::ThreeStatePhi(Real dt) {
    Mat3 Phi_clk{{1, dt, dt * dt / 2}, {0, 1, dt}, {0, 0, 1}};
    return Phi_clk;
  }

  State ClockDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    return PropagateImpl(x0, t0, tf, u, nullptr, nullptr);
  }

  State ClockDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) {
    return PropagateImpl(x0, t0, tf, u, nullptr, stm);
  }

  State ClockDynamics::PropagateWithRelativity(const State& x0, Real t0, Real tf,
                                               const ClockRelativityContext& relativity,
                                               const State* u, MatXd* stm) {
    return PropagateImpl(x0, t0, tf, u, &relativity, stm);
  }

  State ClockDynamics::PropagateImpl(const State& x0, Real t0, Real tf, const State* u,
                                     const ClockRelativityContext* relativity, MatXd* stm) {
    LUPNT_CHECK(u == nullptr, "ClockDynamics does not support control inputs", "ClockDynamics");

    MatXd Phi_clk, Q_clk;
    VecX noise(x0.size());
    if (x0.size() == 2) {
      Phi_clk = TwoStatePhi(tf - t0).cast<double>();
      if (model_ != ClockModel::UNDEFINED && add_noise_) {
        Q_clk = TwoStateNoise(model_, tf - t0, bias_unit_).cast<double>();
        noise = SampleMvNormal(Vec2d::Zero(), Q_clk, 1, rng_.get()).transpose();
      }
    } else if (x0.size() == 3) {
      Phi_clk = ThreeStatePhi(tf - t0).cast<double>();
      if (model_ != ClockModel::UNDEFINED && add_noise_) {
        Q_clk = ThreeStateNoise(model_, tf - t0, bias_unit_).cast<double>();
        noise = SampleMvNormal(Vec3d::Zero(), Q_clk, 1, rng_.get()).transpose();
      }
    } else {
      LUPNT_CHECK(false, fmt::format("Invalid clock state size {}", x0.size()), "ClockDynamics");
    }

    VecX xf = Phi_clk * x0;
    if (relativity != nullptr) {
      Real rate0 = RelativisticRateCorrection(relativity->rv0_centered.head<3>(),
                                              relativity->rv0_centered.tail<3>(), relativity->GM,
                                              relativity->c, relativity->reference_rate_offset);
      Real ratef = RelativisticRateCorrection(relativity->rvf_centered.head<3>(),
                                              relativity->rvf_centered.tail<3>(), relativity->GM,
                                              relativity->c, relativity->reference_rate_offset);
      Real dt = tf - t0;
      Real bias_rel_s = 0.5 * (rate0 + ratef) * dt;
      xf(0) += SecondsToBiasUnits(bias_rel_s, bias_unit_);
    }
    if (model_ != ClockModel::UNDEFINED && add_noise_) xf += noise;
    State x_state(x0, xf);
    x_state.SetUnits(GetStateUnits(x0.size(), bias_unit_));
    if (stm != nullptr) *stm = Phi_clk.cast<double>();
    return x_state;
  }

  REGISTER_FACTORY_CLASS(ClockDynamics, ClockDynamics)

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<ClockDynamics, Config&>::Creator>&
  AssetFactory<ClockDynamics, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<ClockDynamics, Config&>;

}  // namespace lupnt
