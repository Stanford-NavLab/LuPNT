/**
 * @file DynamicsNumerical.cpp
 * @author Stanford NAV LAB
 * @brief List of Numerical Dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/dynamics/numerical_orbit_dynamics.h"

#include <algorithm>
#include <cctype>
#include <limits>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"
#include "lupnt/environment/body.h"
#include "lupnt/environment/forces.h"
#include "lupnt/environment/solar_system.h"
#include "lupnt/interfaces/kernels.h"
#include "lupnt/interfaces/yaml.h"

namespace lupnt {

  using BodyFactory = AssetFactory<Body, YAML::Node&>;

  namespace {
    /// @brief Indirect (frame-origin) acceleration of a third body.
    ///
    /// A frame-relative acceleration is the difference between the acceleration of the
    /// spacecraft and that of the frame's centre. When an attracting body is *not* that
    /// centre it pulls on the centre too, so the direct term must be corrected by
    /// `-GM s / |s|^3`, with `s` the body's position in the integration frame. This is the
    /// second half of the classical third-body formula that `AccelerationPointMass`
    /// already applies; the spherical-harmonic path needs it just as much, and omitting it
    /// there costs `GM / |s|^2` -- for Earth seen from a Moon-centred frame,
    /// 2.7e-3 m/s^2, eight orders of magnitude above the J2 term the harmonics add.
    ///
    /// Returns zero when the body *is* the frame centre (`s = 0`), which is the common
    /// case of a central body carrying its own gravity field.
    Vec3 IndirectTerm(const Vec3& s, Real GM) {
      if (s.norm() <= EPS) return Vec3::Zero();
      return -GM * s / pow(s.norm(), 3);
    }

    std::string LowerUnitName(std::string value) {
      std::transform(value.begin(), value.end(), value.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return value;
    }

    double ParseLengthUnit(const YAML::Node& node) {
      if (!node) return METER;
      if (node.IsScalar()) {
        std::string unit = LowerUnitName(node.as<std::string>());
        if (unit == "m" || unit == "meter" || unit == "meters") return METER;
        if (unit == "km" || unit == "kilometer" || unit == "kilometers") return KILOMETER;
      }
      return node.as<double>();
    }

    double ParseTimeUnit(const YAML::Node& node) {
      if (!node) return SECOND;
      if (node.IsScalar()) {
        std::string unit = LowerUnitName(node.as<std::string>());
        if (unit == "s" || unit == "sec" || unit == "second" || unit == "seconds") return SECOND;
      }
      return node.as<double>();
    }

    double ParseMassUnit(const YAML::Node& node) {
      if (!node) return KILOGRAM;
      if (node.IsScalar()) {
        std::string unit = LowerUnitName(node.as<std::string>());
        if (unit == "kg" || unit == "kilogram" || unit == "kilograms") return KILOGRAM;
      }
      return node.as<double>();
    }

    UnitSystem ParseUnitSystem(const YAML::Node& node) {
      if (!node) return SI_UNITS;
      if (node.IsScalar()) {
        std::string unit = LowerUnitName(node.as<std::string>());
        if (unit == "si" || unit == "m_s_kg" || unit == "mks") return SI_UNITS;
        if (unit == "km_s_kg" || unit == "kms" || unit == "kilometers") return KM_S_KG_UNITS;
      }
      return UnitSystem{ParseLengthUnit(node["length"]), ParseTimeUnit(node["time"]),
                        ParseMassUnit(node["mass"])};
    }

    Vec3 PositionFromSI(const Vec3& r, const UnitSystem& units) { return r / units.length; }
    Vec3 PositionToSI(const Vec3& r, const UnitSystem& units) { return r * units.length; }
    Vec3 VelocityToSI(const Vec3& v, const UnitSystem& units) {
      return v * (units.length / units.time);
    }
    Vec3 AccelerationFromSI(const Vec3& a, const UnitSystem& units) {
      return a * (units.time * units.time / units.length);
    }
    Vec6 StateToSI(const Vec6& rv, const UnitSystem& units) {
      Vec6 out;
      out.head(3) = PositionToSI(rv.head(3), units);
      out.tail(3) = VelocityToSI(rv.tail(3), units);
      return out;
    }
  }  // namespace

  NumericalDynamics::NumericalDynamics() { SetIntegrator(default_integrator); }

  VecX NumericalDynamics::ComputeRates(Real t, const State& x) const { return odefunc_(t, x); }

  NumericalDynamics::NumericalDynamics(Config& config) {
    SetIntegrator(default_integrator);
    if (config["dt"]) SetTimeStep(config["dt"].as<Real>());
    if (config["integrator"]) {
      IntegratorType integ
          = enum_cast<IntegratorType>(config["integrator"].as<std::string>()).value();
      SetIntegrator(integ);
    }
    // Optional adaptive-integrator tolerances / iteration cap. Applied only if at least one
    // key is present, so the integrator's built-in defaults are otherwise preserved.
    if (config["abstol"] || config["reltol"] || config["max_iter"]) {
      IntegratorParams defaults;
      int max_iter = config["max_iter"].as<int>(defaults.max_iter);
      double abstol = config["abstol"].as<double>(defaults.abstol);
      double reltol = config["reltol"].as<double>(defaults.reltol);
      SetIntegratorParams(IntegratorParams(max_iter, abstol, reltol));
    }
    // if (config["params"]) SetParams(config["params"].as<DynamicsParam>());
  }

  void NumericalDynamics::SetTimeStep(Real dt) { dt_ = dt; };
  Real NumericalDynamics::GetTimeStep() const { return dt_; };
  void NumericalDynamics::SetODE(ODE odefunc) { odefunc_ = odefunc; };

  State NumericalDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    if (abs(tf - t0) < EPS) return x0;
    (void)u;
    return State(x0, integrator_->Propagate(odefunc_, t0, tf, x0, dt_));
  }
  // State NumericalOrbitDynamics::Propagate(const State &x0, Real t0, Real tf, const State *u,
  //                                         MatXd *stm) {
  //   if (abs(tf - t0) < EPS) return x0;
  //   (void)u;
  //   return State(x0, integrator_->Propagate(odefunc_, t0, tf, x0, dt_, stm));
  // }
  void NumericalDynamics::SetIntegrator(IntegratorType integ) {
    if (integ == IntegratorType::RK4) {
      integrator_ = MakeUnique<RK4>();
    } else if (integ == IntegratorType::RK8) {
      integrator_ = MakeUnique<RK8>();
    } else if (integ == IntegratorType::RKF45) {
      integrator_ = MakeUnique<RKF45>();
    } else if (integ == IntegratorType::PD45) {
      integrator_ = MakeUnique<PD45>();
    } else {
      LUPNT_CHECK(false, "Invalid Integrator Type", "NumericalOrbitDynamics");
    }
  }

  void NumericalDynamics::SetIntegratorParams(IntegratorParams params) {
    integrator_->SetParams(params);
  }

  // New overloads with termination info (no STM)
  State NumericalDynamics::PropagateEx(const State& x0, Real t0, Real tf, TerminationInfo* info) {
    if (std::abs((tf - t0).val()) < EPS) {
      if (info) {
        info->terminated = false;
        info->reason = TerminationReason::ReachedTf;
        info->t_stop = t0;
        info->steps = 0;
      }
      return x0;
    }
    IntegratorResult res = integrator_->PropagateEx(odefunc_, t0, tf, x0, dt_);
    if (info) {
      info->terminated = (res.reason != TerminationReason::ReachedTf);
      info->reason = res.reason;
      info->t_stop = res.t;
      info->steps = res.steps;
    }
    return State(x0, res.x);
  }

  MatX NumericalDynamics::PropagateEx(const State& x0, Real t0, const VecX& tf,
                                      TerminationInfo* info) {
    MatX xf = MatX::Zero(tf.size(), x0.size());
    Ptr<ProgressBar> pbar = nullptr;
    if (print_progress_) {
      pbar = Logger::GetProgressBar(tf.size(), "Propagating", "Propagation");
      pbar->SetDescription("Propagating");
    }
    xf.row(0) = PropagateEx(x0, t0, tf(0), info);
    for (int i = 1; i < tf.size(); i++) {
      Real t0_i = tf(i - 1);
      Real tf_i = tf(i);
      VecX x0_i = xf.row(i - 1);
      VecX xf_i = PropagateEx(x0_i, t0_i, tf_i, info);
      xf.row(i) = xf_i;
      if (info && info->terminated) {
        if (print_progress_) pbar->Finish();
        return xf.topRows(i + 1);  // return only the completed part
      }
      if (print_progress_) pbar->Update(i);
    }
    if (print_progress_) pbar->Finish();
    return xf;
  }

  // New overloads with termination info + STM
  State NumericalDynamics::PropagateExStm(const State& x0, Real t0, Real tf, MatXd* stm,
                                          TerminationInfo* info) {
    if (std::abs((tf - t0).val()) < EPS) {
      if (info) {
        info->terminated = false;
        info->reason = TerminationReason::ReachedTf;
        info->t_stop = t0;
        info->steps = 0;
      }
      return x0;
    }
    IntegratorResult res = integrator_->PropagateEx(odefunc_, t0, tf, x0, dt_, stm);
    if (info) {
      info->terminated = (res.reason != TerminationReason::ReachedTf);
      info->reason = res.reason;
      info->t_stop = res.t;
      info->steps = res.steps;
    }
    return State(x0, res.x);
  }

  // ****************************************************************************
  // CartesianTwoBodyDynamics
  // ****************************************************************************

  CartesianTwoBodyDynamics::CartesianTwoBodyDynamics(Real GM) : NumericalDynamics(), GM_(GM) {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
  };

  VecX CartesianTwoBodyDynamics::ComputeRates(Real t, const State& x) const {
    (void)t;
    Vec6 rv_dot;
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    rv_dot.head(3) = v;
    rv_dot.tail(3) = -GM_ * r / pow(r_norm, 3);

    return rv_dot;
  }

  // ****************************************************************************
  // JToCartTwoBodyDynamics
  // ****************************************************************************

  JToCartTwoBodyDynamics::JToCartTwoBodyDynamics(Real GM, Real J2, Real R_body, Frame frame,
                                                 Frame body_fixed_frame)
      : NumericalDynamics(),
        GM_(GM),
        J2_(J2),
        R_body_(R_body),
        frame_(frame),
        body_fixed_frame_(body_fixed_frame) {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
  };

  VecX JToCartTwoBodyDynamics::ComputeRates(Real t, const State& x) const {
    Vec6 rv_dot(6);
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    rv_dot.head(3) = v;
    rv_dot.tail(3) = -GM_ * r / pow(r_norm, 3);

    if (abs(J2_) <= EPS) return rv_dot;

    LUPNT_CHECK(frame_ != Frame::UNDEFINED, "Frame not set", "JToCartTwoBodyDynamics");
    LUPNT_CHECK(body_fixed_frame_ != Frame::UNDEFINED, "Body-fixed frame not set",
                "JToCartTwoBodyDynamics");
    LUPNT_CHECK(GetFrameCenter(frame_) == GetFrameCenter(body_fixed_frame_),
                "J2 dynamics requires integration and body-fixed frames with the same origin",
                "JToCartTwoBodyDynamics");

    Real t_tdb = t + GetLupntEpoch();
    auto [R_frame_to_body_fixed, translation]
        = GetFrameRotationTranslation(t_tdb, frame_, body_fixed_frame_);
    (void)translation;

    Vec3 r_bf = R_frame_to_body_fixed * r;
    Vec3 a_J2_bf;
    Real aux1 = -3.0 / 2.0 * GM_ * J2_ * pow(R_body_, 2.0) / pow(r_norm, 5.0);
    Real aux2 = 5.0 * pow(r_bf(2) / r_norm, 2.0);
    a_J2_bf(0) = aux1 * (1.0 - aux2) * r_bf(0);
    a_J2_bf(1) = aux1 * (1.0 - aux2) * r_bf(1);
    a_J2_bf(2) = aux1 * (3.0 - aux2) * r_bf(2);

    rv_dot.tail(3) += R_frame_to_body_fixed.transpose() * a_J2_bf;

    return rv_dot;
  }

  // ****************************************************************************
  // J2KeplerianDynamics
  // ****************************************************************************

  J2KeplerianDynamics::J2KeplerianDynamics(Real GM, Real J2, Real R_body)
      : NumericalDynamics(), GM_(GM), J2_(J2), R_body_(R_body) {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
  };

  VecX J2KeplerianDynamics::ComputeRates(Real t, const State& x) const {
    (void)t;
    Real p = x(0) * (1.0 - x(1) * x(1));
    Real n = sqrt(GM_ / pow(x(0), 3.0));
    Real eta = sqrt(1.0 - x(1) * x(1));

    Vec6 coe_dot(6);
    coe_dot(0) = 0.0;
    coe_dot(1) = 0.0;
    coe_dot(2) = 0.0;
    coe_dot(3) = -3.0 / 2.0 * J2_ * pow(R_body_ / p, 2.0) * n * cos(x(2));
    coe_dot(4) = 3.0 / 4.0 * J2_ * pow(R_body_ / p, 2.0) * n * (5.0 * pow(cos(x(2)), 2.0) - 1.0);
    coe_dot(5)
        = n + 3.0 / 4.0 * J2_ * pow(R_body_ / p, 2.0) * n * eta * (3.0 * pow(cos(x(2)), 2.0) - 1.0);
    return coe_dot;
  }

  // ****************************************************************************
  // MoonMeanDynamics
  // ****************************************************************************

  MoonMeanDynamics::MoonMeanDynamics() : NumericalDynamics() {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
  }

  VecX MoonMeanDynamics::ComputeRates(Real t, const State& x) const {
    (void)t;
    Real a = x(0);
    Real e = x(1);
    Real i = x(2);
    Real w = x(4);

    Vec6 coe_dot(6);
    coe_dot(0) = 0;
    coe_dot(1) = (15 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON)) * e
                 * sqrt(1 - pow(e, 2)) * pow(sin(i), 2) * sin(2 * w);
    coe_dot(2) = -(15 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (16 * sqrt(GM_MOON)) * pow(e, 2)
                 / sqrt(1 - pow(e, 2)) * sin(2 * i) * sin(2 * w);
    coe_dot(3) = -(3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                     / (2 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * cos(i)
                 + (3 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0))
                       / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                       * (5 * pow(e, 2) * cos(2 * w) - 3 * pow(e, 2) - 2) * cos(i);
    coe_dot(4) = (3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                     / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * (5 * pow(cos(i), 2) - 1)
                 + (3 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0))
                       / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                       * ((5 * pow(cos(i), 2) - 1 + pow(e, 2))
                          + 5 * (1 - pow(e, 2) - pow(cos(i), 2)) * cos(2 * w));
    coe_dot(5) = sqrt(GM_MOON) / pow(a, 3.0 / 2.0)
                 + (3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                       / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 3.0 / 2.0))
                       * (3 * pow(cos(i), 2) - 1)
                 - (k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON))
                       * ((3 * pow(e, 2) + 7) * (3 * pow(cos(i), 2) - 1)
                          + 15 * (1 + pow(e, 2)) * pow(sin(i), 2) * cos(2 * w));
    return coe_dot;
  }

  // ****************************************************************************
  // NBodyDynamics
  // ****************************************************************************
  void NBodyDynamics::AddBody(const Body& body) {
    for (auto& b : bodies_) {
      LUPNT_CHECK(b.id != body.id, "Body already added", "NBodyDynamics");
    }
    LUPNT_CHECK(body.units == units_, "Body unit system must match NBodyDynamics unit system",
                "NBodyDynamics");
    bodies_.push_back(body);
  }

  void NBodyDynamics::RemoveBody(const Body& body) {
    for (auto it = bodies_.begin(); it != bodies_.end(); ++it) {
      if (it->id == body.id) {
        bodies_.erase(it);
        break;
      }
    }
  }

  void NBodyDynamics::SetSrpCoeff(Real bcoeff) {
    SetParam("bcoeff_srp", bcoeff);
    use_srp_ = true;
  }

  void NBodyDynamics::SetDragCoeff(Real bcoeff) {
    SetParam("bcoeff_drag", bcoeff);
    use_drag_ = true;
  }

  void NBodyDynamics::SetSrpCoefficient(Real CR, Real area, Real mass) {
    Real bcoeff = (CR * area) / mass;
    SetParam("bcoeff_srp", bcoeff);
    use_srp_ = true;
  }

  void NBodyDynamics::SetDragCoefficient(Real CD, Real area, Real mass) {
    Real bcoeff = (CD * area) / mass;
    SetParam("bcoeff_drag", bcoeff);
    use_drag_ = true;
  }

  NBodyDynamics::NBodyDynamics() : NumericalDynamics() {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
    ParamState default_params = ParamState(VecX::Zero(2), {"bcoeff_srp", "bcoeff_drag"});
    this->SetParams(default_params);
  };

  ForceModelSpec ParseForceModelSpec(const Config& force_model) {
    ForceModelSpec spec;
    if (force_model["bodies"]) {
      for (const auto& body_config : force_model["bodies"]) {
        for (const auto& body : body_config) {
          const std::string body_name = body.first.as<std::string>();
          auto id = enum_cast<BodyId>(body_name);
          LUPNT_CHECK(id.has_value(),
                      fmt::format("Invalid body name in force_model: {}", body_name), "ForceModel");
          switch (id.value()) {
            case BodyId::MOON:
              spec.moon_degree = body.second["n"].as<int>(0);
              spec.moon_order = body.second["m"].as<int>(0);
              break;
            case BodyId::EARTH: spec.include_earth = true; break;
            case BodyId::SUN: spec.include_sun = true; break;
            default: break;
          }
        }
      }
    }
    // Accept both spellings so the shared style matches world/agent (`relativity`) and the
    // per-example filter/truth blocks that historically used `use_relativity`.
    if (force_model["relativity"]) spec.relativity = force_model["relativity"].as<bool>();
    if (force_model["use_relativity"]) spec.relativity = force_model["use_relativity"].as<bool>();
    // SRP ballistic coefficient (optional).
    if (force_model["CR"]) {
      spec.has_srp = true;
      spec.srp_cr = force_model["CR"].as<double>();
      spec.srp_area_m2 = force_model["area"].as<double>(spec.srp_area_m2);
      spec.srp_mass_kg = force_model["mass"].as<double>(spec.srp_mass_kg);
    }
    return spec;
  }

  NBodyDynamics::NBodyDynamics(Config& dynamics_config) : NumericalDynamics(dynamics_config) {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
    Logger::Debug("Creating", "NBodyDynamics");
    ParamState default_params = ParamState(VecX::Zero(2), {"bcoeff_srp", "bcoeff_drag"});
    this->SetParams(default_params);

    if (dynamics_config["units"]) {
      units_ = ParseUnitSystem(dynamics_config["units"]);
      Logger::Debug(fmt::format("Setting unit system length={} m, time={} s, mass={} kg",
                                units_.length, units_.time, units_.mass),
                    "NBodyDynamics");
    }

    // SRP and drag
    if (dynamics_config["area"] && dynamics_config["mass"]) {
      Real mass = dynamics_config["mass"].as<Real>();
      Real area = dynamics_config["area"].as<Real>();
      if (dynamics_config["CR"]) {
        Real CR = dynamics_config["CR"].as<Real>();
        SetSrpCoefficient(CR, area, mass);
        Logger::Debug("Setting SRP Coefficient", "NBodyDynamics");
      }
      if (dynamics_config["CD"]) {
        Real CD = dynamics_config["CD"].as<Real>();
        SetDragCoefficient(CD, area, mass);
        Logger::Debug("Setting Drag Coefficient", "NBodyDynamics");
      }
    }

    // Bodies
    if (dynamics_config["bodies"]) {
      const auto& bodies = dynamics_config["bodies"];
      for (const auto& body_config : bodies) {
        for (const auto& body : body_config) {
          std::string body_name = body.first.as<std::string>();
          Logger::Debug(fmt::format("Adding body {}", body_name), "NBodyDynamics");
          LUPNT_CHECK(enum_cast<BodyId>(body_name).has_value(),
                      fmt::format("Invalid body name: {}", body_name), "NBodyDynamics");
          BodyId body_id = enum_cast<BodyId>(body_name).value();
          int n = body.second["n"].as<int>(0);
          int m = body.second["m"].as<int>(0);
          std::string gravity_file = body.second["gravity_file"].as<std::string>("");
          AddBody(Body::CreateBody(body_id, units_, n, m, gravity_file));
          Logger::Debug(fmt::format("Adding body {}", body_name), "NBodyDynamics");
        }
      }
    }

    // Frame
    if (dynamics_config["frame"]) {
      frame_ = enum_cast<Frame>(dynamics_config["frame"].as<std::string>()).value();
      Logger::Debug(fmt::format("Setting frame to {}", enum_name(frame_)), "NBodyDynamics");
    }

    // Autodiff
    if (dynamics_config["autodiff"]) {
      use_ad_ = dynamics_config["autodiff"].as<bool>();
      Logger::Debug(fmt::format("Setting autodiff to {}", use_ad_), "NBodyDynamics");
    }

    if (dynamics_config["relativity"]) {
      use_relativity_ = dynamics_config["relativity"].as<bool>();
      Logger::Debug(fmt::format("Setting relativity to {}", use_relativity_), "NBodyDynamics");
    }
  }

  REGISTER_FACTORY_CLASS(Dynamics, NBodyDynamics)

  void NBodyDynamics::SetUnits(const UnitSystem& units) {
    LUPNT_CHECK(bodies_.empty(), "Set units before adding bodies", "NBodyDynamics");
    units_ = units;
  }

  State NBodyDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) {
    LUPNT_CHECK(use_ad_, "Autodiff not enabled", "NBodyDynamics");
    (void)u;
    return NumericalDynamics::Propagate(x0, t0, tf, u, stm);
  }

  Vec3 NBodyDynamics::RelativisticNBodyAcceleration(Real t_tdb, const Vec3& r,
                                                    const Vec3& v) const {
    // Moyer (2000) Eq. (4-26) is written in Solar-System barycentric coordinates,
    // so gather SSB-referenced states of the spacecraft and the massive bodies.
    // Position differences are frame-independent, but the absolute velocities in
    // the 1/c^2 terms are not, so the barycentric velocities matter.
    PhysicalConstants constants = GetPhysicalConstants(units_);
    Vec6 rv_center = GetBodyPosVel(t_tdb, BodyId::SSB, GetFrameCenter(frame_), frame_, units_);
    Vec3 r_ssb = r + rv_center.head(3);
    Vec3 v_ssb = v + rv_center.tail(3);
    std::vector<Vec3> r_bodies, v_bodies;
    std::vector<Real> mu_bodies;
    r_bodies.reserve(bodies_.size());
    v_bodies.reserve(bodies_.size());
    mu_bodies.reserve(bodies_.size());
    for (const auto& body : bodies_) {
      Vec6 rv_b = GetBodyPosVel(t_tdb, BodyId::SSB, body.id, frame_, units_);
      r_bodies.push_back(rv_b.head(3));
      v_bodies.push_back(rv_b.tail(3));
      mu_bodies.push_back(body.GM);
    }
    return AccelerationRelativisticNBody(r_ssb, v_ssb, r_bodies, v_bodies, mu_bodies, constants.C);
  }

  Vec3 NBodyDynamics::AccelerationSrp(Real t_tdb, const Vec3& r) const {
    PhysicalConstants constants = GetPhysicalConstants(units_);

    // Radiation pressure from the configured irradiance rather than the built-in default.
    // `solar_flux_` and `C` are SI, so the SI pressure is scaled into `units_` here.
    Real p_sun = units_.Pressure(solar_flux_ / C);

    // The Sun's position must share an origin with `r` (the integration frame's centre) for
    // the Sun->spacecraft direction to be correct.
    Vec3 r_sun = GetBodyPos(t_tdb, BodyId::SUN, frame_, units_);

    // Any body other than the Sun can occult it. Shadowing is multiplicative across
    // occulters, and each body's illumination must be evaluated in that body's own frame.
    Real illumination = 1.0;
    for (const auto& body : bodies_) {
      if (body.id == BodyId::SUN) continue;
      Vec3 r_body = GetBodyPos(t_tdb, body.id, frame_, units_);  // zero for the frame centre
      illumination *= Illumination(r - r_body, r_sun - r_body, body.R);
    }

    return illumination * AccelerationSolarRadiation(r, r_sun, GetSrpCoeff(), p_sun, constants.AU);
  }

  VecX NBodyDynamics::ComputeRates(Real t, const State& rv) const {
    Real t_tdb = t + GetLupntEpoch();
    LUPNT_CHECK(frame_ != Frame::UNDEFINED, "Frame not set", "NBodyDynamics");

    // Position, velocity, and acceleration in the configured coherent unit system.
    // w.r.t. to the inertial frame origin
    Vec3 r = rv.head(3);
    Vec3 v = rv.tail(3);
    Vec3 a = Vec3::Zero();

    for (const auto& body : bodies_) {
      if (body.use_gravity_field) {
        auto& grav = body.gravity_field;
        Vec3 r_si = PositionToSI(r, units_);
        Vec3 r_bf = PositionFromSI(ConvertFrame(t_tdb, r_si, frame_, body.fixed_frame), units_);
        Vec3 a_bf;
        if (use_ad_) {
          a_bf = AccelarationGravityField<Real>(r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m);
        } else {
          Vec3d r_bf_d = r_bf.cast<double>();
          a_bf = AccelarationGravityField<double>(r_bf_d, static_cast<double>(grav.GM),
                                                  static_cast<double>(grav.R),
                                                  grav.CS.cast<double>(), grav.n, grav.m);
        }
        auto [R_bf_to_frame, translation]
            = GetFrameRotationTranslation(t_tdb, body.fixed_frame, frame_);
        (void)translation;
        Vec3 r_body = GetBodyPos(t_tdb, body.id, frame_, units_);
        Vec3 ai = R_bf_to_frame * a_bf + IndirectTerm(r_body, grav.GM);
        a += ai;
      } else {
        Vec3 r_body = GetBodyPos(t_tdb, body.id, frame_, units_);
        Vec3 ai = AccelerationPointMass(rv.head(3), r_body, body.GM);
        a += ai;
      }

      // Atmospheric drag
      if (use_drag_ && body.id == BodyId::EARTH) {
        // TODO: Currently only works for Earth
        Real tt = ConvertTime(t_tdb, Time::TDB, Time::TT);
        Real mjd_tt = TimeToMjd(tt);
        MatX3 Rot = NutationMatrix(mjd_tt) * PrecessionMatrix(MJD_J2000_TT, mjd_tt);
        Vec6 rv_si = StateToSI(rv, units_);
        Real bcoeff_si = units_.ToSI(GetDragCoeff().val(), 2, 0, -1);
        Vec3 a_drag_si = AccelerationDrag(mjd_tt, rv_si, Rot, bcoeff_si);
        a += AccelerationFromSI(a_drag_si, units_);
      }
    }

    // Solar radiation pressure: one term for the spacecraft, not one per body.
    if (use_srp_) a += AccelerationSrp(t_tdb, r);

    if (use_relativity_) {
      a += RelativisticNBodyAcceleration(t_tdb, r, v);
    }

    Vec6 rv_dot;
    rv_dot << v, a;
    return rv_dot;
  }

  std::map<std::string, Vec3> NBodyDynamics::ComputeAccelerations(Real t, const State& rv,
                                                                  bool decompose_gravity) const {
    Real t_tdb = t + GetLupntEpoch();
    LUPNT_CHECK(frame_ != Frame::UNDEFINED, "Frame not set", "NBodyDynamics");

    Vec3 r = rv.head(3);
    Vec3 v = rv.tail(3);

    std::map<std::string, Vec3> acc;

    Vec3 a_drag = Vec3::Zero();

    for (const auto& body : bodies_) {
      if (body.use_gravity_field) {
        auto& grav = body.gravity_field;
        Vec3 r_si = PositionToSI(r, units_);
        Vec3 r_bf = PositionFromSI(ConvertFrame(t_tdb, r_si, frame_, body.fixed_frame), units_);
        // The central (point-mass) term corresponds to the C(0,0)=1 coefficient.
        Real r_bf_norm = r_bf.norm();
        Vec3 a_bf_central = -grav.GM * r_bf / pow(r_bf_norm, 3);
        auto [R_bf_to_frame, translation]
            = GetFrameRotationTranslation(t_tdb, body.fixed_frame, frame_);
        (void)translation;
        // The central term carries the indirect (frame-origin) contribution, so
        // `<body>_gravity` matches what the point-mass branch reports for the same body
        // and the `_Jn`/`_Cnm`/`_nonspherical` keys stay purely non-spherical.
        Vec3 r_body = GetBodyPos(t_tdb, body.id, frame_, units_);
        acc[body.name + "_gravity"] = R_bf_to_frame * a_bf_central + IndirectTerm(r_body, grav.GM);

        if (decompose_gravity) {
          // The gravity-field acceleration is linear in the spherical-harmonic
          // coefficients, so the contribution of a single (n, m) harmonic is
          // recovered by evaluating the field with a coefficient matrix in which
          // only that (n, m) pair (C_nm, and S_nm for m>0) is retained.
          for (int m = 0; m <= grav.m; m++) {
            for (int n = std::max(m, 1); n <= grav.n; n++) {
              MatX CS_iso = MatX::Zero(grav.CS.rows(), grav.CS.cols());
              CS_iso(n, m) = grav.CS(n, m);                      // C_nm
              if (m >= 1) CS_iso(m - 1, n) = grav.CS(m - 1, n);  // S_nm
              Vec3 a_bf_nm
                  = AccelarationGravityField<Real>(r_bf, grav.GM, grav.R, CS_iso, grav.n, grav.m);
              std::string key = (m == 0) ? body.name + "_J" + std::to_string(n)
                                         : body.name + "_C" + std::to_string(n) + std::to_string(m);
              acc[key] = R_bf_to_frame * a_bf_nm;
            }
          }
        } else {
          Vec3 a_bf
              = AccelarationGravityField<Real>(r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m);
          acc[body.name + "_nonspherical"] = R_bf_to_frame * (a_bf - a_bf_central);
        }
      } else {
        Vec3 r_body = GetBodyPos(t_tdb, body.id, frame_, units_);
        acc[body.name + "_gravity"] = AccelerationPointMass(r, r_body, body.GM);
      }

      // Atmospheric drag
      if (use_drag_ && body.id == BodyId::EARTH) {
        Real tt = ConvertTime(t_tdb, Time::TDB, Time::TT);
        Real mjd_tt = TimeToMjd(tt);
        MatX3 Rot = NutationMatrix(mjd_tt) * PrecessionMatrix(MJD_J2000_TT, mjd_tt);
        Vec6 rv_si = StateToSI(rv, units_);
        Real bcoeff_si = units_.ToSI(GetDragCoeff().val(), 2, 0, -1);
        Vec3 a_drag_si = AccelerationDrag(mjd_tt, rv_si, Rot, bcoeff_si);
        a_drag += AccelerationFromSI(a_drag_si, units_);
      }
    }

    // Solar radiation pressure: one term for the spacecraft, not one per body.
    if (use_srp_) acc["srp"] = AccelerationSrp(t_tdb, r);
    if (use_drag_) acc["drag"] = a_drag;

    if (use_relativity_) {
      acc["relativity"] = RelativisticNBodyAcceleration(t_tdb, r, v);
    }

    return acc;
  }
};  // namespace lupnt
