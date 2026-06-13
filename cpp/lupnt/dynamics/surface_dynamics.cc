/**
 * @file dynamics_surface.cc
 * @author Stanford NAV Lab
 * @brief  Dynamics of objects on Surface (Rover, GroundStation)
 * @version 0.1
 * @date 2024-11-26
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/dynamics/surface_dynamics.h"

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/definitions.h"
#include "lupnt/interfaces/yaml.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  // StaticDynamics
  StaticDynamics::StaticDynamics() : Dynamics() {}

  StaticDynamics::StaticDynamics(Config& config) : Dynamics(config) {}

  State StaticDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    (void)t0;
    (void)tf;
    (void)u;
    return x0;
  }

  REGISTER_FACTORY_CLASS(Dynamics, StaticDynamics)

  // Surface2DDynamics
  SurfaceDynamics2D::SurfaceDynamics2D() : Dynamics() {}

  SurfaceDynamics2D::SurfaceDynamics2D(Config& config) : Dynamics(config) {
    if (config["seed"]) {
      seed_ = config["seed"].as<int>();
      rng_ = MakePtr<std::mt19937>(seed_);
    }
    Q_.setZero();
    if (config["sigma_r"]) {
      Real sigma_r = config["sigma_r"].as<Real>();
      Q_(0, 0) = sigma_r * sigma_r;
      Q_(1, 1) = sigma_r * sigma_r;
      add_noise_ = true;
    }
    if (config["sigma_theta"]) {
      Real sigma_theta = RAD * config["sigma_theta"].as<Real>();
      Q_(2, 2) = sigma_theta * sigma_theta;
      add_noise_ = true;
    }
  }

  State SurfaceDynamics2D::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    LUPNT_CHECK(u, "Control is required", "Surface2DDynamics");
    LUPNT_CHECK(u->size() == 2, "Control must be a 2D vector", "Surface2DDynamics");

    Real dt = tf - t0;
    State xf = x0;
    xf(0) += dt * (*u)(0) * cos(x0(2));
    xf(1) += dt * (*u)(0) * sin(x0(2));
    xf(2) += dt * (*u)(1);

    if (add_noise_) xf += SampleMvNormal(Vec3d::Zero(), Q_ * dt * dt, 1, rng_.get()).transpose();

    return xf;
  }

  REGISTER_FACTORY_CLASS(Dynamics, SurfaceDynamics2D)
}  // namespace lupnt
