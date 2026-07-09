#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/filters/filter.h"
#include "lupnt/states/state.h"

namespace lupnt {

  /// @brief Unified measurement *data*: a predicted/observed observable vector and its
  /// noise covariance at a given epoch.
  ///
  /// This is the return type of every `Measurement` and `ErrorStateMeasurement` model's
  /// `Compute`. `value` stacks the observables selected by the
  /// measurement's config (e.g. pseudorange/Doppler rows for GNSS, range/range-rate for a
  /// crosslink), and `covariance` is the corresponding diagonal (by convention) measurement
  /// noise `R`, sized `value.size() x value.size()`.
  struct MeasData {
    Real timestamp = 0.0;  ///< epoch [s]
    VecXd value;           ///< stacked observable vector `z`
    MatXd covariance;      ///< measurement noise covariance `R` [n_z x n_z]
  };

  /// @brief One inertial-measurement-unit sample (gyro + accelerometer), used as an *input*
  /// to strapdown mechanization in `Predict`, not as a measurement-update model.
  ///
  /// Kept as a plain data holder (it does not participate in the `Measurement` model
  /// hierarchy): the mechanization needs the full nav state and attitude, so no model
  /// methods live here.
  struct ImuMeasurement {
    Real timestamp = 0.0;   ///< sample epoch [s]
    Vec3 w = Vec3::Zero();  ///< [rad/s] angular velocity
    Vec3 a = Vec3::Zero();  ///< [m/s^2] acceleration
  };

  /// @brief Abstract base class for a measurement *model* driven by a filter `State` vector.
  ///
  /// A concrete `Measurement` owns its own measurement-model math and configuration
  /// (settings + state-index mapping, typically a nested `Config`). It maps an estimation
  /// `State` to a predicted observable vector `z = h(x)`, optionally with its Jacobian
  /// `H = dh/dx` (`Compute`, with an optional `H` output), and can wrap itself into a
  /// `FilterMeasurementFunction` for direct use with the `Filter` family (`CreateFunction`).
  ///
  /// This is the base for state-vector measurements (GNSS, inter-satellite crosslinks,
  /// ground-station range/range-rate, etc.). Error-state (INS) measurements, whose Jacobian
  /// is taken with respect to a nav *error* state and depends on attitude, use the parallel
  /// `ErrorStateMeasurement` base instead.
  class Measurement {
  public:
    virtual ~Measurement() = default;

    /// @brief Predict the measurement `z = h(x)`, its noise covariance `R`, and (if
    /// `H != nullptr`) the measurement Jacobian `H = dh/dx` sized `n_z x x.size()`.
    /// @param x      Estimation state.
    /// @param[out] H If non-null, filled with the measurement Jacobian; if null, only the
    ///               value and covariance are computed.
    /// @return `MeasData` with `value = z` and `covariance = R`.
    virtual MeasData Compute(const State& x, MatXd* H = nullptr) const = 0;

    /// @brief Polymorphic deep copy so `CreateFunction`'s closure can own an independent
    /// copy of this model (avoids dangling references / slicing).
    virtual Ptr<Measurement> Clone() const = 0;

    /// @brief Wrap this model into a `FilterMeasurementFunction` `(x, H, R) -> z` suitable
    /// for `Filter::SetMeasurementFunction`.
    ///
    /// The default implementation captures a `Clone()` of this model and forwards to
    /// `Compute(x, H)`, copying the returned `MeasData::covariance` into `R`.
    /// Subclasses may override for specialized covariance/Jacobian handling.
    virtual FilterMeasurementFunction CreateFunction() const {
      Ptr<Measurement> self = Clone();
      return [self](const State& x, MatXd* H, MatXd* R) -> VecXd {
        MeasData md = self->Compute(x, H);
        if (R != nullptr) *R = md.covariance;
        return md.value;
      };
    }
  };

  /// @brief CRTP helper that supplies `Clone()` for a concrete `Measurement` subclass.
  ///
  /// Usage: `class GnssMeasurement : public MeasurementClone<GnssMeasurement> { ... };`
  /// (still virtual through the `Measurement` base). Requires the derived type to be
  /// copy-constructible.
  template <typename Derived> class MeasurementClone : public Measurement {
  public:
    Ptr<Measurement> Clone() const override {
      return MakePtr<Derived>(static_cast<const Derived&>(*this));
    }
  };

  /// @brief Nominal navigation state an `ErrorStateMeasurement` model linearizes about.
  ///
  /// Error-state (INS) filters keep the nominal state as structured members and estimate a
  /// small error state `[dr, dv, dtheta, db_a, db_g, d(clock_bias), d(clock_drift)]`. Because
  /// the measurement Jacobian depends on the attitude `R_b2n` (which does not live naturally
  /// in a flat `State` vector), the nominal state is passed to the model through this context
  /// rather than as a `State`.
  struct NavErrorContext {
    Vec3d r = Vec3d::Zero();          ///< position, Moon-fixed frame [m]
    Vec3d v = Vec3d::Zero();          ///< velocity, Moon-fixed frame [m/s]
    Mat3d R_b2n = Mat3d::Identity();  ///< body-to-nav attitude
    double clock_bias_s = 0.0;        ///< receiver clock bias [s]
    double clock_drift_sps = 0.0;     ///< receiver clock drift [s/s]
    int error_state_size = 0;         ///< dimension of the error state (columns of `H`)
  };

  /// @brief Abstract base class for an error-state (INS) measurement *model*.
  ///
  /// Concrete subclasses own their measurement-model math and configuration (settings +
  /// error-state index mapping) and carry the measured observable(s) they compare against.
  /// Unlike `Measurement`, the Jacobian is taken with respect to the nav *error* state and
  /// depends on the nominal attitude, so the model consumes a `NavErrorContext` instead of a
  /// `State`, and there is no `FilterMeasurementFunction`/`Filter` integration (the hosting
  /// application drives the Joseph-form error-state update directly).
  class ErrorStateMeasurement {
  public:
    virtual ~ErrorStateMeasurement() = default;

    /// @brief Predict the measurement `z = h(nominal)`, its noise covariance `R`, and (if
    /// `H != nullptr`) the Jacobian `H = dh/d(error state)` sized
    /// `n_z x nominal.error_state_size`.
    /// @param nominal Nominal nav state to linearize about.
    /// @param[out] H  If non-null, filled with the error-state measurement Jacobian.
    /// @return `MeasData` with `value = z_pred` and `covariance = R`.
    virtual MeasData Compute(const NavErrorContext& nominal, MatXd* H = nullptr) const = 0;
  };

}  // namespace lupnt
