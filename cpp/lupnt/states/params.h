/**
 * @file state.h
 * @author Stanford NAV LAB
 * @brief  Paramter state class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once
#include "lupnt/states/state.h"

namespace lupnt {

  /// @brief Labeled vector of "extra" dynamics parameters (e.g. SRP/drag
  /// coefficients, integer ambiguities) carried alongside a primary `State`,
  /// always expressed in frame MOON_CI.
  ///
  /// `ParamState` is the parameter container threaded through
  /// `Dynamics::GetParams`/`Dynamics::SetParams` and
  /// `Dynamics::PropagateWithParams`, and is concatenated across sub-states by
  /// `JointState::Add`/`JointState::GetParams` so the EKF/UKF/batch filters
  /// can estimate, consider, or fix these parameters alongside the orbit/clock
  /// state.
  class ParamState : public State {
  public:
    /// @brief Construct an empty parameter state (type `Params`, frame MOON_CI).
    ParamState() : State() {
      SetType(TYPE);
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct a zero-initialized parameter state of length `n`
    /// (type `Params`, frame MOON_CI), with placeholder names/units.
    /// @param n  Number of parameters.
    ParamState(int n) : State(n) {
      SetType(TYPE);
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct a parameter state from values and names, inferring
    /// units from the known parameter names (see GetParamUnits).
    ///
    /// Used e.g. by `NumericalOrbitDynamics::GetParams` to build the default
    /// `{"bcoeff_srp", "bcoeff_drag"}` parameter state.
    ///
    /// @param x      Parameter values.
    /// @param names  Parameter names (same length as `x`), e.g.
    ///               `"bcoeff_srp"`, `"bcoeff_drag"`, `"IntegerAmbiguity"`.
    ParamState(const VecX& x, std::vector<std::string> names) : State(x) {
      SetType(TYPE);
      SetNames(names);
      SetUnits(GetParamUnits(names));
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct a parameter state from values, names, and explicit units.
    ///
    /// Used by `Dynamics::PropagateWithParams` to rebuild a `ParamState` after
    /// propagation, preserving the original parameter names/units.
    ///
    /// @param x      Parameter values.
    /// @param names  Parameter names (same length as `x`).
    /// @param units  Parameter unit strings (same length as `x`).
    ParamState(const VecX& x, std::vector<std::string> names, std::vector<std::string> units)
        : State(x) {
      SetType(TYPE);
      SetNames(names);
      SetUnits(units);
      SetFrame(Frame::MOON_CI);
    }

    static constexpr StateType TYPE = "Params";

  private:
    /// @brief Infer per-parameter unit strings from their names, for use by
    /// the `ParamState(x, names)` constructor when no explicit units are
    /// given.
    ///
    /// Recognizes `"bcoeff_drag"`/`"bcoeff_srp"` (drag/SRP area-to-mass
    /// coefficients, [m^2/kg]); all other names (including
    /// `"IntegerAmbiguity"`) default to dimensionless `"-"`.
    ///
    /// @param names  Parameter names.
    /// @return       Unit string per name (or `{"-"}` if `names` is empty).
    std::vector<std::string> GetParamUnits(const std::vector<std::string>& names) {
      if (names.empty()) {
        return {"-"};
      }
      std::vector<std::string> units;
      for (const auto& name : names) {
        if (name == "bcoeff_drag" || name == "bcoeff_srp") {
          units.push_back("m^2/kg");
        } else if (name == "IntegerAmbiguity") {
          units.push_back("-");
        } else {
          units.push_back("-");
        }
      }
      return units;
    }
  };

}  // namespace lupnt
