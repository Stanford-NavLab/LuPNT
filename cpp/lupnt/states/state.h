#pragma once
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/error.h"

namespace lupnt {

  using StateType = std::string_view;

  // State *********************************************************************

  /// @brief Generic labeled state vector: a `VecX` augmented with a type tag,
  /// per-element names/units, and a reference `Frame`.
  ///
  /// This is the common base class for all concrete state representations used
  /// throughout LuPNT (e.g. `Cart6`, `ClassicalOE`, `ClockState3`,
  /// `AttitudeState`). Dynamics models (`Dynamics::Propagate`,
  /// `Dynamics::PropagateWithParams`), frame conversions (`ConvertFrame`,
  /// `ConvertState`), filters (`Filter`, `JointState`), and measurement models
  /// pass and return `State` objects so that the numerical vector always
  /// carries along its type/frame/units metadata for bookkeeping, logging
  /// (`ToString`), and runtime type checks (`LUPNT_CHECK(x.GetType() == ...)`).
  class State : public VecX {
  private:
    std::string name_ = "State";
    StateType type_ = UNDEFINED;
    std::vector<std::string> names_;
    std::vector<std::string> units_;
    Frame frame_;

    /// @brief Resize the per-element `names_`/`units_` vectors to length `n`
    /// and fill them with placeholder values ("0", "1", ... and "none").
    ///
    /// Called by the size-only constructors (`State(int)`, `State(int, int)`)
    /// so that a freshly-sized, type-less `State` still has well-formed
    /// names/units arrays before a derived class (e.g. `Cart6::Init`)
    /// overwrites them with its real labels.
    void SetNamesAndUnits(int n);

  public:
    /// @brief Construct an empty (zero-length) state.
    State();

    /// @brief Construct a zero-initialized state of length `n`, with
    /// placeholder names/units (see SetNamesAndUnits).
    ///
    /// @param n  Number of elements in the state vector.
    State(int n);

    /// @brief Construct a zero-initialized column-vector state with `rows`
    /// elements; `cols` must be 1.
    ///
    /// @param rows  Number of elements (vector length).
    /// @param cols  Must be 1 (states are always column vectors); checked via
    ///              LUPNT_CHECK.
    State(int rows, int cols);

    /// @brief Construct a state from a raw vector together with explicit
    /// per-element `names` and `units`.
    ///
    /// Used e.g. by `Dynamics::Propagate`/`filter_utils` to wrap a propagated
    /// numerical vector back into a `State` while preserving the original
    /// state's element labels and units.
    ///
    /// @param vec    Underlying numerical values.
    /// @param names  Per-element variable names (same length as `vec`).
    /// @param units  Per-element unit strings (same length as `vec`).
    State(const VecX& vec, const std::vector<std::string>& names,
          const std::vector<std::string>& units);

    /// @brief Construct a state by copying the type/names/units/frame metadata
    /// from `state` but using the numerical values from `vec`.
    ///
    /// Used when a derived-state type needs to rebuild itself from a plain
    /// Eigen vector (e.g. after an arithmetic operation) while keeping the
    /// original state's labeling.
    ///
    /// @param state  Source state to copy type/names/units/frame from.
    /// @param vec    New numerical values (must match `state`'s size).
    template <typename T, int N> State(const State& state, const Vector<T, N>& vec);

    /// @brief Copy constructor: copies the numerical vector plus
    /// type/names/units/frame metadata from `other`.
    State(const State& other);

    /// @brief Construct a state directly from any Eigen column- or row-vector
    /// expression, assigning placeholder names/units via SetNamesAndUnits.
    ///
    /// Enables implicit/explicit construction of `State` (and derived types)
    /// from Eigen vector expressions (e.g. `Vec6`, `head()`/`tail()` slices)
    /// throughout the dynamics and conversion code.
    ///
    /// @param x  Eigen vector (or row-vector) expression to copy values from.
    template <typename Derived, typename = std::enable_if_t<(Derived::ColsAtCompileTime == 1
                                                             || Derived::RowsAtCompileTime == 1)>>
    State(const Eigen::DenseBase<Derived>& x) {
      if (x.cols() == 1) {
        VecX::operator=(x);
        SetNamesAndUnits(x.rows());
      } else if (x.rows() == 1) {
        VecX::operator=(x.transpose());
        SetNamesAndUnits(x.cols());
      } else {
        static_assert(Derived::ColsAtCompileTime == 1 || Derived::RowsAtCompileTime == 1,
                      "State can only be constructed from vectors");
      }
    }

    /// @brief Set the reference frame (e.g. Frame::GCRF, Frame::MOON_CI) this
    /// state's coordinates are expressed in. Read by `ConvertFrame` /
    /// `ConvertState` to determine the source frame for a conversion.
    void SetFrame(Frame frame);

    /// @brief Set the state-type tag (e.g. `Cart6::TYPE`, `ClassicalOE::TYPE`),
    /// used by `LUPNT_CHECK(x.GetType() == ...)` runtime checks throughout the
    /// dynamics/conversion code to verify a `State` holds the expected
    /// representation.
    void SetType(StateType type);

    /// @brief Set the human-readable name of this state (e.g. for logging).
    void SetName(const std::string& name);

    /// @brief Set the per-element variable names (e.g. {"r_x","r_y",...}),
    /// used by `ToString`/formatters and by `JointState` when concatenating
    /// sub-states into a combined labeled vector.
    void SetNames(const std::vector<std::string>& names);

    /// @brief Set the per-element unit strings (e.g. {"m","m","m","m/s",...}),
    /// used by `ToString`/formatters and by `JointState` when concatenating
    /// sub-states into a combined labeled vector.
    void SetUnits(const std::vector<std::string>& units);

    /// @brief Get the reference frame this state's coordinates are expressed in.
    /// @return Reference frame (e.g. Frame::GCRF, Frame::MOON_CI).
    Frame GetFrame() const;

    /// @brief Get the state-type tag (e.g. `Cart6::TYPE`).
    /// @return State type identifier string.
    StateType GetType() const;

    /// @brief Get the human-readable name of this state.
    std::string GetName() const;

    /// @brief Get the per-element variable names.
    /// @return Names, one per element of this state vector.
    std::vector<std::string> GetNames() const;

    /// @brief Get the per-element unit strings.
    /// @return Units, one per element of this state vector.
    std::vector<std::string> GetUnits() const;

    /// @brief Copy-assign the numerical vector plus type/names/units/frame
    /// metadata from `other`.
    State& operator=(const State& other);

    /// @brief Assign this state's numerical values from any Eigen column- or
    /// row-vector expression, leaving existing type/names/units/frame metadata
    /// unchanged.
    ///
    /// @param x  Eigen vector (or row-vector) expression to copy values from
    ///           (must match this state's size).
    template <typename Derived, typename = std::enable_if_t<(Derived::ColsAtCompileTime == 1
                                                             || Derived::RowsAtCompileTime == 1)>>
    State& operator=(const Eigen::DenseBase<Derived>& x) {
      if (x.cols() == 1) {
        VecX::operator=(x);
      } else if (x.rows() == 1) {
        VecX::operator=(x.transpose());
      } else {
        static_assert(Derived::ColsAtCompileTime == 1 || Derived::RowsAtCompileTime == 1,
                      "State can only be assigned from vectors");
      }
      return *this;
    }

    static constexpr StateType TYPE = UNDEFINED;

    /// @brief Format this state as `"<Type>(<values>, <Frame>)"` for logging
    /// and debugging (e.g. via `fmt::format` / the `fmt::formatter<State>`
    /// specialization below).
    /// @return Human-readable string representation of this state.
    std::string ToString() const;
  };

  template <typename T, int N> State::State(const State& state, const Vector<T, N>& vec)
      : VecX(vec) {
    type_ = state.type_;
    names_ = state.names_;
    units_ = state.units_;
    frame_ = state.frame_;
  }

  // Cart6 *********************************************************************
  // - $r_x$ [m] (position x-component)
  // - $r_y$ [m] (position y-component)
  // - $r_z$ [m] (position z-component)
  // - $v_x$ [m/s] (velocity x-component)
  // - $v_y$ [m/s] (velocity y-component)
  // - $v_z$ [m/s] (velocity z-component)
  class Cart6 : public State {
  private:
    void Init() {
      SetType(Cart6::TYPE);
      SetNames({"r_x", "r_y", "r_z", "v_x", "v_y", "v_z"});
      SetUnits({"m", "m", "m", "m/s", "m/s", "m/s"});
    }

  public:
    static constexpr StateType TYPE = "Cart6";

    /// @brief Reinterpret a generic `State` as a `Cart6` (Cartesian
    /// position+velocity), checking its size and type tag.
    ///
    /// Used wherever dynamics/conversion code receives a `State` and needs to
    /// access it as `r()`/`v()` (e.g. force-model accelerations, frame
    /// conversions via `ConvertState`).
    ///
    /// @param x  6-element state of type `Cart6::TYPE`.
    Cart6(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 6, "State size must be 6", "Cart6");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be Cart6", "Cart6");
      Init();
    }

    /// @brief Construct a `Cart6` from separate position and velocity vectors.
    /// @param r      Position [m] in `frame`.
    /// @param v      Velocity [m/s] in `frame`.
    /// @param frame  Reference frame (default MOON_CI).
    Cart6(const Vec3& r, const Vec3& v, const Frame frame = Frame::MOON_CI) : State() {
      Init();
      SetFrame(frame);
      resize(6);
      head(3) = r;
      tail(3) = v;
    }

    /// @brief Construct a `Cart6` from a stacked 6-element [position; velocity]
    /// vector.
    /// @param x      [r_x, r_y, r_z, v_x, v_y, v_z] in [m, m, m, m/s, m/s, m/s].
    /// @param frame  Reference frame (default MOON_CI).
    Cart6(const Vec6& x = Vec6::Zero(), const Frame frame = Frame::MOON_CI) : State(x) {
      Init();
      SetFrame(frame);
    }

    /// @brief Mutable view of the position sub-vector.
    /// @return Position [m] in this state's frame.
    Eigen::Ref<Vec3> r() { return head(3); }

    /// @brief Mutable view of the velocity sub-vector.
    /// @return Velocity [m/s] in this state's frame.
    Eigen::Ref<Vec3> v() { return tail(3); }

    /// @brief Position sub-vector. @return Position [m] in this state's frame.
    Vec3 r() const { return head(3); }

    /// @brief Velocity sub-vector. @return Velocity [m/s] in this state's frame.
    Vec3 v() const { return tail(3); }
  };

  // AttitudeState ***************************************************************
  // - $q_0$ [-] (scalar part)
  // - $q_1$ [-] (x-component)
  // - $q_2$ [-] (y-component)
  // - $q_3$ [-] (z-component)
  // - $w_x$ [rad/s] (angular velocity x-component)
  // - $w_y$ [rad/s] (angular velocity y-component)
  // - $w_z$ [rad/s] (angular velocity z-component)
  class Attitude : public State {
  private:
    void Init() {
      resize(7);
      head(4) = Vec4(1.0, 0.0, 0.0, 0.0);
      tail(3) = Vec3(0.0, 0.0, 0.0);
      SetType("AttitudeState");
      SetNames({"q0", "q1", "q2", "q3", "wx", "wy", "wz"});
      SetUnits({"-", "-", "-", "-", "rad/s", "rad/s", "rad/s"});
    }

  public:
    static constexpr StateType TYPE = "AttitudeState";

    /// @brief Construct an identity-attitude state: zero rotation
    /// (q = [1,0,0,0]) and zero angular velocity.
    Attitude() : State() { Init(); }

    /// @brief Construct a zero-initialized attitude state of the given shape;
    /// `cols` must be 1 (see `State(int, int)`).
    Attitude(int rows, int cols) : State(rows, cols) { Init(); }

    /// @brief Reinterpret a generic `State` as an `Attitude` state, checking
    /// its size and type tag.
    /// @param x  7-element state of type `Attitude::TYPE` (quaternion + angular velocity).
    Attitude(const State& x) : Attitude() {
      LUPNT_CHECK(x.size() == 7, "State size must be 7", "AttitudeState");
      LUPNT_CHECK(x.GetType() == "AttitudeState", "State type must be AttitudeState",
                  "AttitudeState");
      head(7) = x.head(7);
      SetFrame(x.GetFrame());
    }

    /// @brief Construct from a stacked 7-element [quaternion; angular velocity] vector.
    /// @param qw  [q0, q1, q2, q3, wx, wy, wz] in [-, -, -, -, rad/s, rad/s, rad/s].
    Attitude(const Vec7& qw) : Attitude() { head(7) = qw; }

    /// @brief Construct from a separate attitude quaternion and angular velocity.
    /// @param q  Attitude quaternion [-] (scalar-first: q0, q1, q2, q3).
    /// @param w  Angular velocity [rad/s].
    /// @param f  Reference frame the quaternion/angular velocity are expressed in.
    Attitude(const Vec4& q, const Vec3& w, Frame f) : Attitude() {
      head(4) = q;
      tail(3) = w;
      SetFrame(f);
    }

    /// @brief Mutable view of the attitude quaternion. @return Quaternion [-] (q0,q1,q2,q3).
    Eigen::Ref<Vec4> q() { return head(4); }

    /// @brief Mutable view of the angular velocity. @return Angular velocity [rad/s].
    Eigen::Ref<Vec3> w() { return tail(3); }

    /// @brief Attitude quaternion. @return Quaternion [-] (q0,q1,q2,q3).
    Vec4 q() const { return head(4); }

    /// @brief Angular velocity. @return Angular velocity [rad/s].
    Vec3 w() const { return tail(3); }
  };

  // RelCart6 *******************************************************************
  // - $\delta r_x$ [m] (position x-component)
  // - $\delta r_y$ [m] (position y-component)
  // - $\delta r_z$ [m] (position z-component)
  // - $\delta v_x$ [m/s] (velocity x-component)
  // - $\delta v_y$ [m/s] (velocity y-component)
  // - $\delta v_z$ [m/s] (velocity z-component)
  class RelCart6 : public State {
  private:
    void Init() {
      resize(6);
      SetType(Cart6::TYPE);
      SetNames({"dr_x", "dr_y", "dr_z", "dv_x", "dv_y", "dv_z"});
      SetUnits({"m", "m", "m", "m/s", "m/s", "m/s"});
    }

  public:
    static constexpr StateType TYPE = "RelCart6";

    /// @brief Reinterpret a generic `State` as a `RelCart6` (relative
    /// Cartesian position+velocity, e.g. chief-to-deputy in RTN/synodic
    /// frames), checking its size and type tag. Used by `StateConverter` when
    /// converting between `Cart6` and relative-orbit representations
    /// (`InertialToSynodic` and ROE conversions).
    /// @param x  6-element state of type `RelCart6::TYPE`.
    RelCart6(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 6, "State size must be 6", "RelCart6");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be RelCart6", "RelCart6");
      Init();
    }

    /// @brief Construct a `RelCart6` from separate relative position and velocity vectors.
    /// @param dr     Relative position [m] in `frame`.
    /// @param dv     Relative velocity [m/s] in `frame`.
    /// @param frame  Reference frame (default MOON_CI).
    RelCart6(const Vec3& dr, const Vec3& dv, const Frame frame = Frame::MOON_CI) : State() {
      Init();
      SetFrame(frame);
      resize(6);
      head(3) = dr;
      tail(3) = dv;
    }

    /// @brief Construct a `RelCart6` from a stacked 6-element [relative position; relative
    /// velocity] vector.
    /// @param x      [dr_x, dr_y, dr_z, dv_x, dv_y, dv_z] in [m, m, m, m/s, m/s, m/s].
    /// @param frame  Reference frame (default MOON_CI).
    RelCart6(const Vec6& x = Vec6::Zero(), const Frame frame = Frame::MOON_CI) : State(x) {
      Init();
      SetFrame(frame);
    }

    /// @brief Mutable view of the relative position. @return Relative position [m].
    Eigen::Ref<Vec3> dr() { return head(3); }

    /// @brief Mutable view of the relative velocity. @return Relative velocity [m/s].
    Eigen::Ref<Vec3> dv() { return tail(3); }

    /// @brief Relative position. @return Relative position [m].
    Vec3 dr() const { return head(3); }

    /// @brief Relative velocity. @return Relative velocity [m/s].
    Vec3 dv() const { return tail(3); }
  };

  // SurfaceState2D *************************************************************
  // - $x$ [m] (x-coordinate)
  // - $y$ [m] (y-coordinate)
  // - $\theta$ [rad] (angle)
  class SurfaceState2D : public State {
  private:
    void Init() {
      resize(3);
      SetType(SurfaceState2D::TYPE);
      SetNames({"x", "y", "theta"});
      SetUnits({"m", "m", "rad"});
    }

  public:
    /// @brief Reinterpret a generic `State` as a `SurfaceState2D` (2D
    /// position + heading of a surface rover/lander), checking its size and
    /// type tag.
    /// @param x  3-element state of type `SurfaceState2D::TYPE`.
    SurfaceState2D(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 3, "State size must be 3", "Surface2DState");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be Surface2D", "Surface2DState");
      Init();
    }

    /// @brief Construct a `SurfaceState2D` from a stacked [x, y, theta] vector.
    /// @param x      [x, y, theta] in [m, m, rad].
    /// @param frame  Reference frame (default MOON_CI).
    SurfaceState2D(const Vec3& x = Vec3::Zero(), const Frame frame = Frame::MOON_CI) : State(x) {
      Init();
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "Surface2D";
  };

  // Cart3 **********************************************************************
  // - $r_x$ [m] (position x-component)
  // - $r_y$ [m] (position y-component)
  // - $r_z$ [m] (position z-component)
  class Cart3 : public State {
  private:
    void Init() {
      resize(3);
      SetType(Cart3::TYPE);
      SetNames({"r_x", "r_y", "r_z"});
      SetUnits({"m", "m", "m"});
    }

  public:
    /// @brief Reinterpret a generic `State` as a `Cart3` (3D Cartesian
    /// position only), checking its size and type tag.
    /// @param x  3-element state of type `Cart3::TYPE`.
    Cart3(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 3, "State size must be 3", "Cart3");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be Cart3", "Cart3");
      Init();
    }

    /// @brief Construct a `Cart3` from a 3-element position vector.
    /// @param x      Position [m].
    /// @param frame  Reference frame (default MOON_CI).
    Cart3(const Vec3& x = Vec3::Zero(), const Frame frame = Frame::MOON_CI) : State(x) {
      Init();
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "Cart3";
  };

  // Cart2 **********************************************************************
  // - $x$ [m] (x-coordinate)
  // - $y$ [m] (y-coordinate)
  class Cart2 : public State {
  private:
    void Init() {
      resize(2);
      SetType(Cart2::TYPE);
      SetNames({"x", "y"});
      SetUnits({"m", "m"});
    }

  public:
    /// @brief Reinterpret a generic `State` as a `Cart2` (2D Cartesian
    /// position only), checking its size and type tag.
    /// @param x  2-element state of type `Cart2::TYPE`.
    Cart2(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 2, "State size must be 2", "CartTo");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be CartTo", "CartTo");
      Init();
    }

    /// @brief Construct a `Cart2` from a 2-element position vector.
    /// @param x      Position [m].
    /// @param frame  Reference frame (default MOON_CI).
    Cart2(const Vec2& x = Vec2::Zero(), const Frame frame = Frame::MOON_CI) : State(x) {
      Init();
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "CartTo";

    /// @brief Mutable view of the position. @return Position [m].
    Eigen::Ref<Vec2> r() { return head(2); }

    /// @brief Position. @return Position [m].
    Vec2 r() const { return head(2); }
  };

  // AzElRange ******************************************************************
  // - $az$ [rad] (azimuth)
  // - $el$ [rad] (elevation)
  // - $range$ [m] (range)
  class AzElRange : public State {
  private:
    void Init() {
      resize(3);
      SetType(AzElRange::TYPE);
      SetNames({"az", "el", "range"});
      SetUnits({"rad", "rad", "m"});
    }

  public:
    /// @brief Reinterpret a generic `State` as an `AzElRange` (ground-station
    /// look-angle measurement state), checking its size and type tag. Used by
    /// the coordinate conversions `EastNorthUpToAzElRange`/`CartToAzElRange`
    /// and their inverses when building az/el/range measurement states.
    /// @param x  3-element state of type `AzElRange::TYPE`.
    AzElRange(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 3, "State size must be 3", "AzElRange");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be AzElRange", "AzElRange");
      Init();
    }

    /// @brief Construct an `AzElRange` from a stacked [az, el, range] vector.
    /// @param x      [az, el, range] in [rad, rad, m].
    /// @param frame  Reference frame (default MOON_CI).
    AzElRange(const Vec3& x = Vec3::Zero(), const Frame frame = Frame::MOON_CI) : State() {
      Init();
      head(3) = x;
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "AzElRange";
  };  // namespace lupnt

  // LatLonAlt ******************************************************************
  // - $lat$ [rad] (latitude)
  // - $lon$ [rad] (longitude)
  // - $alt$ [m] (altitude)
  class LatLonAlt : public State {
  private:
    void Init() {
      resize(3);
      SetType(LatLonAlt::TYPE);
      SetNames({"lat", "lon", "alt"});
      SetUnits({"rad", "rad", "m"});
    }

  public:
    /// @brief Reinterpret a generic `State` as a `LatLonAlt` (geodetic
    /// latitude/longitude/altitude state), checking its size and type tag.
    /// Used by `CartToLatLonAlt`/`LatLonAltToCart` when converting body-fixed
    /// Cartesian positions to/from geodetic coordinates (e.g. for ground
    /// stations and surface assets).
    /// @param x  3-element state of type `LatLonAlt::TYPE`.
    LatLonAlt(const State& x) : State() {
      LUPNT_CHECK(x.size() == 3, "State size must be 3", "LatLonAlt");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be LatLonAlt", "LatLonAlt");
      Init();
      head(3) = x.head(3);
    }

    /// @brief Construct a `LatLonAlt` from a stacked [lat, lon, alt] vector.
    /// @param x      [lat, lon, alt] in [rad, rad, m].
    /// @param frame  Reference frame (default MOON_CI).
    LatLonAlt(const Vec3& x = Vec3::Zero(), const Frame frame = Frame::MOON_CI) : State() {
      Init();
      head(3) = x;
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "LatLonAlt";
  };

  // ClassicalOE *****************************************************************
  // Singular at $e \in {0,1}$, and $i \in {0, \pi}$
  // - $a$ [m] (semi-major axis)
  // - $e$ [-] (eccentricity)
  // - $i$ [rad] (inclination)
  // - $\Omega$ [rad] (right ascension of the ascending node)
  // - $\omega$ [rad] (argument of periapsis)
  // - $M$ [rad] (mean anomaly)
  class ClassicalOE : public State {
  private:
    void Init() {
      SetType(ClassicalOE::TYPE);
      SetNames({"a", "e", "i", "Omega", "w", "M"});
      SetUnits({"m", "-", "rad", "rad", "rad", "rad"});
    }

  public:
    /// @brief Reinterpret a generic `State` as `ClassicalOE` (classical
    /// Keplerian orbital elements), checking its size and type tag. Used by
    /// `KeplerianDynamics<ClassicalOE>::Propagate` and the
    /// mean/osculating conversions (`MeanToOsculating`/`OsculatingToMean`)
    /// when operating on an orbit-element state.
    /// @param x  6-element state of type `ClassicalOE::TYPE`.
    ClassicalOE(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 6, "State size must be 6", "ClassicalOE");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be ClassicalOE", "ClassicalOE");
      Init();
    }

    /// @brief Construct a `ClassicalOE` from a stacked 6-element element vector.
    /// @param x      [a, e, i, Omega, w, M] in [m, -, rad, rad, rad, rad].
    /// @param frame  Reference frame the elements are defined in (default MOON_CI).
    ClassicalOE(const Vec6& x = Vec6::Zero(), const Frame frame = Frame::MOON_CI) : State(x) {
      Init();
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "ClassicalOE";

    /// @brief Mutable reference to the semi-major axis. @return Semi-major axis [m].
    Real& a() { return this->coeffRef(0); }
    /// @brief Mutable reference to the eccentricity. @return Eccentricity [-].
    Real& e() { return this->coeffRef(1); }
    /// @brief Mutable reference to the inclination. @return Inclination [rad].
    Real& i() { return this->coeffRef(2); }
    /// @brief Mutable reference to the RAAN. @return Right ascension of the ascending node [rad].
    Real& Omega() { return this->coeffRef(3); }
    /// @brief Mutable reference to the argument of periapsis. @return Argument of periapsis [rad].
    Real& w() { return this->coeffRef(4); }
    /// @brief Mutable reference to the mean anomaly. @return Mean anomaly [rad].
    Real& M() { return this->coeffRef(5); }

    /// @brief Semi-major axis. @return Semi-major axis [m].
    Real a() const { return this->coeff(0); }
    /// @brief Eccentricity. @return Eccentricity [-].
    Real e() const { return this->coeff(1); }
    /// @brief Inclination. @return Inclination [rad].
    Real i() const { return this->coeff(2); }
    /// @brief Right ascension of the ascending node. @return RAAN [rad].
    Real Omega() const { return this->coeff(3); }
    /// @brief Argument of periapsis. @return Argument of periapsis [rad].
    Real w() const { return this->coeff(4); }
    /// @brief Mean anomaly. @return Mean anomaly [rad].
    Real M() const { return this->coeff(5); }
  };

  // QuasiNonsingularOE **********************************************************
  // - $a$ [m] (semi-major axis)
  // - $u$ [-] (mean argument of latitude)
  // - $e_x$ [-] (eccentricity x-component)
  // - $e_y$ [-] (eccentricity y-component)
  // - $i$ [rad] (inclination)
  // - $\Omega$ [rad] (right ascension of the ascending node)
  // - $\omega$ [rad] (argument of periapsis)
  class QuasiNonsingularOE : public State {
  private:
    void Init() {
      SetType(QuasiNonsingularOE::TYPE);
      SetNames({"a", "u", "ex", "ey", "i", "Omega"});
      SetUnits({"m", "-", "-", "-", "rad", "rad"});
    }

  public:
    /// @brief Reinterpret a generic `State` as `QuasiNonsingularOE`
    /// (quasi-nonsingular orbital elements), checking its size and type tag.
    /// Used by `KeplerianDynamics<QuasiNonsingularOE>::Propagate` and the
    /// element-set conversions in `StateConverter`. Avoids the singularities
    /// of `ClassicalOE` at e=0 and i=0.
    /// @param x  6-element state of type `QuasiNonsingularOE::TYPE`.
    QuasiNonsingularOE(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 6, "State size must be 6", "QuasiNonsingularOE");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be QuasiNonsingularOE",
                  "QuasiNonsingularOE");
      Init();
    }

    /// @brief Construct a `QuasiNonsingularOE` from a stacked 6-element element vector.
    /// @param x  [a, u, ex, ey, i, Omega] in [m, -, -, -, rad, rad].
    QuasiNonsingularOE(const Vec6& x = Vec6::Zero()) : State(x) { Init(); }

    /// @brief Construct a `QuasiNonsingularOE` from a stacked element vector and reference frame.
    /// @param x      [a, u, ex, ey, i, Omega] in [m, -, -, -, rad, rad].
    /// @param frame  Reference frame the elements are defined in.
    QuasiNonsingularOE(const Vec6& x, const Frame frame) : QuasiNonsingularOE(x) {
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "QuasiNonsingularOE";

    /// @brief Mutable reference to the semi-major axis. @return Semi-major axis [m].
    Real& a() { return this->coeffRef(0); }
    /// @brief Mutable reference to the mean argument of latitude. @return Mean argument of latitude
    /// [-].
    Real& u() { return this->coeffRef(1); }
    /// @brief Mutable reference to the eccentricity x-component. @return Eccentricity x-component
    /// [-].
    Real& ex() { return this->coeffRef(2); }
    /// @brief Mutable reference to the eccentricity y-component. @return Eccentricity y-component
    /// [-].
    Real& ey() { return this->coeffRef(3); }
    /// @brief Mutable reference to the inclination. @return Inclination [rad].
    Real& i() { return this->coeffRef(4); }
    /// @brief Mutable reference to the RAAN. @return Right ascension of the ascending node [rad].
    Real& Omega() { return this->coeffRef(5); }

    /// @brief Semi-major axis. @return Semi-major axis [m].
    Real a() const { return this->coeff(0); }
    /// @brief Mean argument of latitude. @return Mean argument of latitude [-].
    Real u() const { return this->coeff(1); }
    /// @brief Eccentricity x-component. @return Eccentricity x-component [-].
    Real ex() const { return this->coeff(2); }
    /// @brief Eccentricity y-component. @return Eccentricity y-component [-].
    Real ey() const { return this->coeff(3); }
    /// @brief Inclination. @return Inclination [rad].
    Real i() const { return this->coeff(4); }
    /// @brief Right ascension of the ascending node. @return RAAN [rad].
    Real Omega() const { return this->coeff(5); }

    /// @brief State-type tag. @return `"QuasiNonsingularOE"`.
    StateType GetType() const { return "QuasiNonsingularOE"; }
  };

  // DelaunayOE ******************************************************************
  // - $l$ [rad] (mean longitude)
  // - $g$ [rad] (longitude of periapsis)
  // - $h$ [rad] (longitude of ascending node)
  // - $L$ [rad]
  // - $G$ [rad]
  // - $H$ [rad]
  class DelaunayOE : public State {
  private:
    void Init() {
      SetType(DelaunayOE::TYPE);
      SetNames({"l", "g", "h", "L", "G", "H"});
      SetUnits({"rad", "rad", "rad", "rad", "rad", "rad"});
    }

  public:
    /// @brief Construct a `DelaunayOE` (Delaunay orbital elements) from a
    /// stacked 6-element element vector. Used by element-set conversions in
    /// `StateConverter` for canonical Hamiltonian orbit representations.
    /// @param x  [l, g, h, L, G, H], all in [rad].
    DelaunayOE(const Vec6& x = Vec6::Zero()) : State(x) { Init(); }

    /// @brief Construct a `DelaunayOE` from a stacked element vector and reference frame.
    /// @param x      [l, g, h, L, G, H], all in [rad].
    /// @param frame  Reference frame the elements are defined in.
    DelaunayOE(const Vec6& x, const Frame frame) : DelaunayOE(x) {
      Init();
      SetFrame(frame);
    }
    static constexpr StateType TYPE = "DelaunayOE";

    /// @brief Mutable reference to the mean longitude. @return Mean longitude [rad].
    Real& l() { return this->coeffRef(0); }
    /// @brief Mutable reference to the longitude of periapsis. @return Longitude of periapsis
    /// [rad].
    Real& g() { return this->coeffRef(1); }
    /// @brief Mutable reference to the longitude of ascending node. @return Longitude of ascending
    /// node [rad].
    Real& h() { return this->coeffRef(2); }
    /// @brief Mutable reference to the L Delaunay momentum. @return L [rad].
    Real& L() { return this->coeffRef(3); }
    /// @brief Mutable reference to the G Delaunay momentum. @return G [rad].
    Real& G() { return this->coeffRef(4); }
    /// @brief Mutable reference to the H Delaunay momentum. @return H [rad].
    Real& H() { return this->coeffRef(5); }

    /// @brief Mean longitude. @return Mean longitude [rad].
    Real l() const { return this->coeff(0); }
    /// @brief Longitude of periapsis. @return Longitude of periapsis [rad].
    Real g() const { return this->coeff(1); }
    /// @brief Longitude of ascending node. @return Longitude of ascending node [rad].
    Real h() const { return this->coeff(2); }
    /// @brief L Delaunay momentum. @return L [rad].
    Real L() const { return this->coeff(3); }
    /// @brief G Delaunay momentum. @return G [rad].
    Real G() const { return this->coeff(4); }
    /// @brief H Delaunay momentum. @return H [rad].
    Real H() const { return this->coeff(5); }
  };

  // EquinoctialOE ***************************************************************
  // - $a$ [m] (semi-major axis)
  // - $h$ [-]
  // - $k$ [-]
  // - $p$ [-]
  // - $q$ [-]
  // - $\lambda$ [rad]
  class EquinoctialOE : public State {
  public:
    /// @brief Construct an `EquinoctialOE` (equinoctial orbital elements,
    /// singularity-free at e=0 and i=0) from a stacked 6-element element
    /// vector, in frame MOON_CI. Used by
    /// `KeplerianDynamics<EquinoctialOE>::Propagate` and element-set
    /// conversions in `StateConverter`.
    /// @param x  [a, h, k, p, q, lon] in [m, -, -, -, -, rad].
    EquinoctialOE(const Vec6& x = Vec6::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"a", "h", "k", "p", "q", "lon"});
      SetUnits({"m", "-", "-", "-", "-", "rad"});
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct an `EquinoctialOE` from a stacked element vector and reference frame.
    /// @param x      [a, h, k, p, q, lon] in [m, -, -, -, -, rad].
    /// @param frame  Reference frame the elements are defined in.
    EquinoctialOE(const Vec6& x, const Frame frame) : EquinoctialOE(x) { SetFrame(frame); }
    static constexpr StateType TYPE = "EquinoctialOE";

    /// @brief Mutable reference to the semi-major axis. @return Semi-major axis [m].
    Real& a() { return this->coeffRef(0); }
    /// @brief Mutable reference to the h equinoctial element. @return h [-].
    Real& h() { return this->coeffRef(1); }
    /// @brief Mutable reference to the k equinoctial element. @return k [-].
    Real& k() { return this->coeffRef(2); }
    /// @brief Mutable reference to the p equinoctial element. @return p [-].
    Real& p() { return this->coeffRef(3); }
    /// @brief Mutable reference to the q equinoctial element. @return q [-].
    Real& q() { return this->coeffRef(4); }
    /// @brief Mutable reference to the mean longitude. @return Mean longitude [rad].
    Real& lon() { return this->coeffRef(5); }

    /// @brief Semi-major axis. @return Semi-major axis [m].
    Real a() const { return this->coeff(0); }
    /// @brief h equinoctial element. @return h [-].
    Real h() const { return this->coeff(1); }
    /// @brief k equinoctial element. @return k [-].
    Real k() const { return this->coeff(2); }
    /// @brief p equinoctial element. @return p [-].
    Real p() const { return this->coeff(3); }
    /// @brief q equinoctial element. @return q [-].
    Real q() const { return this->coeff(4); }
    /// @brief Mean longitude. @return Mean longitude [rad].
    Real lon() const { return this->coeff(5); }
  };

  // SingularROE *****************************************************************
  // - $a\delta a$ [m] (semi-major axis)
  // - $a\delta M$ [m] (mean anomaly)
  // - $a\delta e$ [m] (eccentricity)
  // - $a\delta \omega$ [m] (argument of periapsis)
  // - $a\delta i$ [m] (inclination)
  // - $a\delta \Omega$ [m] (right ascension of the ascending node)
  class SingularROE : public State {
  public:
    /// @brief Construct a `SingularROE` (singular relative orbital elements,
    /// scaled by semi-major axis) from a stacked 6-element vector, in frame
    /// MOON_CI. Used by `RoeGeometricMappingDynamics` and `StateConverter` to
    /// represent the relative motion of a deputy with respect to a chief in
    /// terms of differential orbital elements.
    /// @param x  [a*da, a*dM, a*de, a*dw, a*di, a*dOmega], all in [m].
    SingularROE(const Vec6& x = Vec6::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"ada", "adM", "ade", "adw", "adi", "adOmega"});
      SetUnits({"m", "m", "m", "m", "m", "m"});
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct a `SingularROE` from a stacked element vector and reference frame.
    /// @param x      [a*da, a*dM, a*de, a*dw, a*di, a*dOmega], all in [m].
    /// @param frame  Reference frame the elements are defined in.
    SingularROE(const Vec6& x, const Frame frame) : SingularROE(x) { SetFrame(frame); }
    static constexpr StateType TYPE = "SingularROE";

    /// @brief Mutable reference to a*(delta semi-major axis). @return a*da [m].
    Real& ada() { return this->coeffRef(0); }
    /// @brief Mutable reference to a*(delta mean anomaly). @return a*dM [m].
    Real& adM() { return this->coeffRef(1); }
    /// @brief Mutable reference to a*(delta eccentricity). @return a*de [m].
    Real& ade() { return this->coeffRef(2); }
    /// @brief Mutable reference to a*(delta argument of periapsis). @return a*dw [m].
    Real& adw() { return this->coeffRef(3); }
    /// @brief Mutable reference to a*(delta inclination). @return a*di [m].
    Real& adi() { return this->coeffRef(4); }
    /// @brief Mutable reference to a*(delta RAAN). @return a*dOmega [m].
    Real& adOmega() { return this->coeffRef(5); }

    /// @brief a*(delta semi-major axis). @return a*da [m].
    Real ada() const { return this->coeff(0); }
    /// @brief a*(delta mean anomaly). @return a*dM [m].
    Real adM() const { return this->coeff(1); }
    /// @brief a*(delta eccentricity). @return a*de [m].
    Real ade() const { return this->coeff(2); }
    /// @brief a*(delta argument of periapsis). @return a*dw [m].
    Real adw() const { return this->coeff(3); }
    /// @brief a*(delta inclination). @return a*di [m].
    Real adi() const { return this->coeff(4); }
    /// @brief a*(delta RAAN). @return a*dOmega [m].
    Real adOmega() const { return this->coeff(5); }
  };

  // QuasiNonsingROE ************************************************************
  // - $a\delta a$ [m] (semi-major axis)
  // - $a\delta l$ [m] (mean longitude)
  // - $a\delta e_x$ [m] (eccentricity x-component)
  // - $a\delta e_y$ [m] (eccentricity y-component)
  // - $a\delta i_x$ [m] (inclination x-component)
  // - $a\delta i_y$ [m] (inclination y-component)
  class QuasiNonsingROE : public State {
  public:
    /// @brief Construct a `QuasiNonsingROE` (quasi-nonsingular relative
    /// orbital elements, scaled by semi-major axis) from a stacked 6-element
    /// vector, in frame MOON_CI. Singularity-free counterpart of
    /// `SingularROE`, used by the same relative-orbit dynamics/conversion
    /// pipelines.
    /// @param x  [a*da, a*dl, a*dex, a*dey, a*dix, a*diy], all in [m].
    QuasiNonsingROE(const Vec6& x = Vec6::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"ada", "adl", "adex", "adey", "adix", "adiy"});
      SetUnits({"m", "m", "m", "m", "m", "m"});
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct a `QuasiNonsingROE` from a stacked element vector and reference frame.
    /// @param x      [a*da, a*dl, a*dex, a*dey, a*dix, a*diy], all in [m].
    /// @param frame  Reference frame the elements are defined in.
    QuasiNonsingROE(const Vec6& x, const Frame frame) : QuasiNonsingROE(x) { SetFrame(frame); }
    static constexpr StateType TYPE = "QuasiNonsingROE";

    /// @brief Mutable reference to a*(delta semi-major axis). @return a*da [m].
    Real& ada() { return this->coeffRef(0); }
    /// @brief Mutable reference to a*(delta mean longitude). @return a*dl [m].
    Real& adl() { return this->coeffRef(1); }
    /// @brief Mutable reference to a*(delta eccentricity x-component). @return a*dex [m].
    Real& adex() { return this->coeffRef(2); }
    /// @brief Mutable reference to a*(delta eccentricity y-component). @return a*dey [m].
    Real& adey() { return this->coeffRef(3); }
    /// @brief Mutable reference to a*(delta inclination x-component). @return a*dix [m].
    Real& adix() { return this->coeffRef(4); }
    /// @brief Mutable reference to a*(delta inclination y-component). @return a*diy [m].
    Real& adiy() { return this->coeffRef(5); }

    /// @brief a*(delta semi-major axis). @return a*da [m].
    Real ada() const { return this->coeff(0); }
    /// @brief a*(delta mean longitude). @return a*dl [m].
    Real adl() const { return this->coeff(1); }
    /// @brief a*(delta eccentricity x-component). @return a*dex [m].
    Real adex() const { return this->coeff(2); }
    /// @brief a*(delta eccentricity y-component). @return a*dey [m].
    Real adey() const { return this->coeff(3); }
    /// @brief a*(delta inclination x-component). @return a*dix [m].
    Real adix() const { return this->coeff(4); }
    /// @brief a*(delta inclination y-component). @return a*diy [m].
    Real adiy() const { return this->coeff(5); }
  };

  // RollPitchYaw ****************************************************************
  // - $roll$ [rad] (roll)
  // - $pitch$ [rad] (pitch)
  // - $yaw$ [rad] (yaw)
  class RollPitchYaw : public State {
  public:
    /// @brief Construct a `RollPitchYaw` attitude state from a stacked
    /// 3-element Euler-angle vector. Used as an alternative attitude
    /// representation to `Attitude`'s quaternion for reporting/conversion.
    /// @param x  [roll, pitch, yaw], all in [rad].
    RollPitchYaw(const Vec3& x = Vec3::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"roll", "pitch", "yaw"});
      SetUnits({"rad", "rad", "rad"});
    }
    static constexpr StateType TYPE = "RollPitchYaw";

    /// @brief Mutable reference to the roll angle. @return Roll [rad].
    Real& roll() { return this->coeffRef(0); }
    /// @brief Mutable reference to the pitch angle. @return Pitch [rad].
    Real& pitch() { return this->coeffRef(1); }
    /// @brief Mutable reference to the yaw angle. @return Yaw [rad].
    Real& yaw() { return this->coeffRef(2); }

    /// @brief Roll angle. @return Roll [rad].
    Real roll() const { return this->coeff(0); }
    /// @brief Pitch angle. @return Pitch [rad].
    Real pitch() const { return this->coeff(1); }
    /// @brief Yaw angle. @return Yaw [rad].
    Real yaw() const { return this->coeff(2); }
  };

  // Quaternion ****************************************************************
  // - $q_0$ [-] (scalar part)
  // - $q_1$ [-] (x-component)
  // - $q_2$ [-] (y-component)
  // - $q_3$ [-] (z-component)
  class Quaternion : public State {
  public:
    /// @brief Construct a `Quaternion` state from a stacked 4-element
    /// scalar-first quaternion vector. Used as a standalone attitude
    /// representation (e.g. for filter/measurement states distinct from the
    /// combined `Attitude` quaternion+angular-velocity state).
    /// @param x  [q0, q1, q2, q3] (scalar-first), all dimensionless [-].
    Quaternion(const Vec4& x = Vec4::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"q_0", "q_1", "q_2", "q_3"});
      SetUnits({"-", "-", "-", "-"});
    }
    static constexpr StateType TYPE = "Quaternion";

    /// @brief Mutable reference to the quaternion scalar part. @return q0 [-].
    Real& q0() { return this->coeffRef(0); }
    /// @brief Mutable reference to the quaternion x-component. @return q1 [-].
    Real& q1() { return this->coeffRef(1); }
    /// @brief Mutable reference to the quaternion y-component. @return q2 [-].
    Real& q2() { return this->coeffRef(2); }
    /// @brief Mutable reference to the quaternion z-component. @return q3 [-].
    Real& q3() { return this->coeffRef(3); }
  };

  // ImuState *******************************************************************
  // - $b\omega_x$ [rad/s] (body frame angular velocity x-component)
  // - $b\omega_y$ [rad/s] (body frame angular velocity y-component)
  // - $b\omega_z$ [rad/s] (body frame angular velocity z-component)
  // - $a_x$ [m/s^2] (body frame acceleration x-component)
  // - $a_y$ [m/s^2] (body frame acceleration y-component)
  class ImuState : public State {
  private:
    void Init() {
      resize(SIZE);
      setZero();
      SetType(TYPE);
      SetNames({"b_w_x", "b_w_y", "b_w_z", "b_a_x", "b_a_y", "b_a_z"});
      SetUnits({"rad/s", "rad/s", "rad/s", "m/s^2", "m/s^2", "m/s^2"});
    }

    /// @brief Check that `x` has the size and type tag expected for an
    /// `ImuState` (6-element, type `ImuState::TYPE`).
    void Check(const State& x) {
      LUPNT_CHECK(x.size() == SIZE, "State size must be 6", "ImuState");
      LUPNT_CHECK(x.GetType() == TYPE, "State type must be ImuState", "ImuState");
    }

  public:
    /// @brief Construct an `ImuState` (IMU bias/measurement state: body-frame
    /// angular velocity + acceleration) from a stacked 6-element vector. Used
    /// by `ImuDynamics::Propagate` as the state propagated for an IMU device.
    /// @param x  [b_wx, b_wy, b_wz, b_ax, b_ay, b_az] in [rad/s, rad/s, rad/s, m/s^2, m/s^2,
    /// m/s^2].
    ImuState(const Vec6& x = Vec6::Zero()) : State() {
      Init();
      head(SIZE) = x;
    }

    /// @brief Reinterpret a generic `State` as an `ImuState`, checking its
    /// size and type tag (see Check).
    /// @param x  6-element state of type `ImuState::TYPE`.
    ImuState(const State& x) : State() {
      Check(x);
      Init();
      head(SIZE) = x.head(SIZE);
    }

    /// @brief Assign from a generic `State`, checking its size and type tag
    /// (see Check) and copying its values.
    ImuState& operator=(const State& x) {
      if (this != &x) {
        Check(x);
        Init();
        head(SIZE) = x.head(SIZE);
      }
      return *this;
    }

    /// @brief Mutable view of the body-frame angular velocity. @return Angular velocity [rad/s].
    Eigen::Ref<Vec3> b_w() { return head(3); }
    /// @brief Mutable view of the body-frame acceleration. @return Acceleration [m/s^2].
    Eigen::Ref<Vec3> b_a() { return tail(3); }
    /// @brief Body-frame angular velocity. @return Angular velocity [rad/s].
    Vec3 b_w() const { return head(3); }
    /// @brief Body-frame acceleration. @return Acceleration [m/s^2].
    Vec3 b_a() const { return tail(3); }

    static constexpr StateType TYPE = "ImuState";
    static constexpr int SIZE = 6;
  };

  // ImuMeasurement *************************************************************
  // - $w_x$ [rad/s] (angular velocity x-component)
  // - $w_y$ [rad/s] (angular velocity y-component)
  // - $w_z$ [rad/s] (angular velocity z-component)
  class AngularVelocity : public State {
  public:
    /// @brief Construct an `AngularVelocity` state from a stacked 3-element
    /// vector. Used as a standalone IMU measurement / output state (e.g.
    /// gyroscope readings) distinct from `ImuState`'s combined bias state.
    /// @param x  [wx, wy, wz], all in [rad/s].
    AngularVelocity(const Vec3& x = Vec3::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"w_x", "w_y", "w_z"});
      SetUnits({"rad/s", "rad/s", "rad/s"});
    }
    static constexpr StateType TYPE = "AngularVelocity";

    /// @brief Mutable reference to the x-component. @return Angular velocity x-component [rad/s].
    Real& wx() { return this->coeffRef(0); }
    /// @brief Mutable reference to the y-component. @return Angular velocity y-component [rad/s].
    Real& wy() { return this->coeffRef(1); }
    /// @brief Mutable reference to the z-component. @return Angular velocity z-component [rad/s].
    Real& wz() { return this->coeffRef(2); }

    /// @brief Angular velocity x-component. @return [rad/s].
    Real wx() const { return this->coeff(0); }
    /// @brief Angular velocity y-component. @return [rad/s].
    Real wy() const { return this->coeff(1); }
    /// @brief Angular velocity z-component. @return [rad/s].
    Real wz() const { return this->coeff(2); }
  };

  class Acceleration : public State {
  public:
    /// @brief Construct an `Acceleration` state from a stacked 3-element
    /// vector. Used as a standalone IMU measurement / output state (e.g.
    /// accelerometer readings) distinct from `ImuState`'s combined bias state.
    /// @param x  [ax, ay, az], all in [m/s^2].
    Acceleration(const Vec3& x = Vec3::Zero()) : State(x) {
      SetType(TYPE);
      SetNames({"a_x", "a_y", "a_z"});
      SetUnits({"m/s^2", "m/s^2", "m/s^2"});
    }
    static constexpr StateType TYPE = "Acceleration";

    /// @brief Mutable reference to the x-component. @return Acceleration x-component [m/s^2].
    Real& ax() { return this->coeffRef(0); }
    /// @brief Mutable reference to the y-component. @return Acceleration y-component [m/s^2].
    Real& ay() { return this->coeffRef(1); }
    /// @brief Mutable reference to the z-component. @return Acceleration z-component [m/s^2].
    Real& az() { return this->coeffRef(2); }

    /// @brief Acceleration x-component. @return [m/s^2].
    Real ax() const { return this->coeff(0); }
    /// @brief Acceleration y-component. @return [m/s^2].
    Real ay() const { return this->coeff(1); }
    /// @brief Acceleration z-component. @return [m/s^2].
    Real az() const { return this->coeff(2); }
  };

  // ClockState ******************************************************************
  class ClockState3 : public State {
  private:
    void Init() {
      resize(3);
      setZero();
      SetType("ClockState");
      SetNames({"b", "d", "dr"});
      SetUnits({"s", "s/s", "s/s^2"});
      SetFrame(Frame::MOON_CI);
    }

  public:
    /// @brief Construct a zero-initialized 3-state clock state
    /// (bias, drift, drift-rate), in frame MOON_CI.
    ClockState3() : State() { Init(); }

    /// @brief Construct a zero-initialized 3-state clock state of the given
    /// shape; `cols` must be 1 (see `State(int, int)`).
    ClockState3(int rows, int cols) : State(rows, cols) { Init(); }

    /// @brief Reinterpret a generic `State` as a `ClockState3`, checking its
    /// size. Used by `Clock::SetState`/`Clock::GetState` and
    /// `ClockDynamics::Propagate` as the receiver/transmitter clock state for
    /// a `Device`.
    /// @param x  3-element state [b, d, dr].
    ClockState3(const State& x) : State() {
      LUPNT_CHECK(x.size() == 3, "State must be of size 2 or 3", "ClockState");
      Init();
      head(3) = x.head(3);
    }

    /// @brief Assign from a generic `State`, checking its size and copying
    /// the first 3 elements.
    ClockState3& operator=(const State& x) {
      if (this != &x) {
        LUPNT_CHECK(x.size() == 3, "State must be of size 2 or 3", "ClockState");
        head(3) = x.head(3);
      }
      return *this;
    }

    /// @brief Mutable reference to the clock bias. @return Clock bias [s].
    Real& b() { return this->coeffRef(0); }
    /// @brief Mutable reference to the clock drift. @return Clock drift [s/s].
    Real& d() { return this->coeffRef(1); }
    /// @brief Mutable reference to the clock drift-rate. @return Clock drift-rate [s/s^2].
    Real& dr() { return this->coeffRef(2); }

    /// @brief Clock bias. @return Clock bias [s].
    Real b() const { return this->coeff(0); }
    /// @brief Clock drift. @return Clock drift [s/s].
    Real d() const { return this->coeff(1); }
    /// @brief Clock drift-rate. @return Clock drift-rate [s/s^2].
    Real dr() const { return this->coeff(2); }
  };

  // ClockState ******************************************************************
  class ClockState2 : public State {
  private:
    void Init() {
      resize(2);
      setZero();
      SetType("ClockState");
      SetNames({"b", "d"});
      SetUnits({"s", "s/s"});
      SetFrame(Frame::MOON_CI);
    }

  public:
    /// @brief Construct a zero-initialized 2-state clock state (bias, drift),
    /// in frame MOON_CI. Used e.g. as a transmitter clock-bias/drift state
    /// carried by a measurement.
    ClockState2() : State() { Init(); }

    /// @brief Construct a zero-initialized 2-state clock state of the given
    /// shape; `cols` must be 1 (see `State(int, int)`).
    ClockState2(int rows, int cols) : State(rows, cols) { Init(); }

    /// @brief Reinterpret a generic `State` as a `ClockState2`, checking its size.
    /// @param x  2-element state [b, d].
    ClockState2(const State& x) : State() {
      LUPNT_CHECK(x.size() == 2, "State must be of size 2 or 3", "ClockState");
      Init();
      head(2) = x.head(2);
    }

    /// @brief Assign from a generic `State`, checking its size and copying
    /// the first 2 elements.
    ClockState2& operator=(const State& x) {
      if (this != &x) {
        LUPNT_CHECK(x.size() == 2, "State must be of size 2 or 3", "ClockState");
        head(2) = x.head(2);
      }
      return *this;
    }

    /// @brief Mutable reference to the clock bias. @return Clock bias [s].
    Real& b() { return this->coeffRef(0); }
    /// @brief Mutable reference to the clock drift. @return Clock drift [s/s].
    Real& d() { return this->coeffRef(1); }

    /// @brief Clock bias. @return Clock bias [s].
    Real b() const { return this->coeff(0); }
    /// @brief Clock drift. @return Clock drift [s/s].
    Real d() const { return this->coeff(1); }
  };

  // JointOrbitClockState ********************************************************
  // - Cartesian orbit state followed by a 2- or 3-state clock model.

  /// @brief Combined Cartesian orbit (position+velocity) and 2- or 3-state
  /// clock state, stacked as [r, v, b, d, (dr)].
  ///
  /// This is the state type propagated by `JointOrbitClockDynamics`, which
  /// integrates a spacecraft's orbital dynamics together with its onboard
  /// clock model in a single state vector (e.g. for GNSS/measurement
  /// pipelines that need a self-consistent orbit + clock-bias estimate).
  class JointOrbitClockState : public State {
  private:
    static constexpr int ORBIT_STATE_SIZE = 6;
    int clock_state_size_ = 3;

    /// @brief Resize/relabel this state for the given clock-state size (2 or
    /// 3) with the provided clock units, setting names to
    /// {"r_x",...,"v_z","b","d",["dr"]} and units to
    /// {"m","m","m","m/s","m/s","m/s", <clock_units>}.
    void Init(int clock_size, const std::vector<std::string>& clock_units) {
      LUPNT_CHECK(clock_size == 2 || clock_size == 3, "Clock state size must be 2 or 3",
                  "JointOrbitClockState");
      clock_state_size_ = clock_size;
      resize(ORBIT_STATE_SIZE + clock_state_size_);
      SetType(JointOrbitClockState::TYPE);
      std::vector<std::string> names = {"r_x", "r_y", "r_z", "v_x", "v_y", "v_z", "b", "d"};
      if (clock_state_size_ == 3) names.push_back("dr");
      SetNames(names);

      std::vector<std::string> units = {"m", "m", "m", "m/s", "m/s", "m/s"};
      units.insert(units.end(), clock_units.begin(), clock_units.end());
      SetUnits(units);
    }

    /// @brief Determine the clock-portion unit strings to use for
    /// `clock_state`: its own units if fully labeled, otherwise the default
    /// {"s","s/s"} or {"s","s/s","s/s^2"} based on its size.
    std::vector<std::string> GetClockUnitsFromState(const State& clock_state) const {
      std::vector<std::string> units = clock_state.GetUnits();
      if (units.size() == static_cast<size_t>(clock_state.size())) return units;
      return clock_state.size() == 2 ? std::vector<std::string>{"s", "s/s"}
                                     : std::vector<std::string>{"s", "s/s", "s/s^2"};
    }

  public:
    static constexpr StateType TYPE = "JointOrbitClock";

    /// @brief Construct a zero-initialized state with a 3-state clock
    /// (bias, drift, drift-rate), in frame MOON_CI.
    JointOrbitClockState() : State() {
      setZero(9);
      Init(3, {"s", "s/s", "s/s^2"});
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Construct a zero-initialized state with a clock of the given
    /// size (2 or 3 states), in frame MOON_CI.
    /// @param clock_size  Number of clock states: 2 (bias, drift) or 3
    ///                     (bias, drift, drift-rate).
    explicit JointOrbitClockState(int clock_size) : State() {
      setZero(ORBIT_STATE_SIZE + clock_size);
      Init(clock_size, clock_size == 2 ? std::vector<std::string>{"s", "s/s"}
                                       : std::vector<std::string>{"s", "s/s", "s/s^2"});
      SetFrame(Frame::MOON_CI);
    }

    /// @brief Reinterpret a generic `State` as a `JointOrbitClockState`,
    /// checking its size (8 or 9) and (re-)labeling it if it isn't already
    /// fully named/unit-tagged. Used by `JointOrbitClockDynamics::Propagate`
    /// to wrap the raw propagated vector.
    /// @param x  8- or 9-element state: [r(3), v(3), b, d, (dr)].
    JointOrbitClockState(const State& x) : State(x) {
      LUPNT_CHECK(x.size() == 8 || x.size() == 9, "State size must be 8 or 9",
                  "JointOrbitClockState");
      clock_state_size_ = x.size() - ORBIT_STATE_SIZE;
      SetType(JointOrbitClockState::TYPE);
      if (x.GetNames().size() != static_cast<size_t>(x.size())
          || x.GetUnits().size() != static_cast<size_t>(x.size())) {
        Init(clock_state_size_, clock_state_size_ == 2
                                    ? std::vector<std::string>{"s", "s/s"}
                                    : std::vector<std::string>{"s", "s/s", "s/s^2"});
        head(x.size()) = x.head(x.size());
        SetFrame(x.GetFrame());
      }
    }

    /// @brief Construct a `JointOrbitClockState` by concatenating a separate
    /// 6-element orbit state and a 2- or 3-element clock state. Used to build
    /// the initial joint state from `Cart6` + `ClockState2`/`ClockState3`
    /// before propagation by `JointOrbitClockDynamics`.
    /// @param orbit_state  6-element Cartesian [r; v] state.
    /// @param clock_state  2- or 3-element clock state [b, d, (dr)].
    JointOrbitClockState(const State& orbit_state, const State& clock_state) : State() {
      LUPNT_CHECK(orbit_state.size() == ORBIT_STATE_SIZE, "Orbit state size must be 6",
                  "JointOrbitClockState");
      LUPNT_CHECK(clock_state.size() == 2 || clock_state.size() == 3,
                  "Clock state size must be 2 or 3", "JointOrbitClockState");
      Init(clock_state.size(), GetClockUnitsFromState(clock_state));
      head(ORBIT_STATE_SIZE) = orbit_state.head(ORBIT_STATE_SIZE);
      tail(clock_state.size()) = clock_state.head(clock_state.size());
      SetFrame(orbit_state.GetFrame());
    }

    /// @brief Assign from a generic `State`, checking its size (8 or 9) and
    /// copying its values and frame.
    JointOrbitClockState& operator=(const State& x) {
      LUPNT_CHECK(x.size() == 8 || x.size() == 9, "State size must be 8 or 9",
                  "JointOrbitClockState");
      if (this != &x) {
        clock_state_size_ = x.size() - ORBIT_STATE_SIZE;
        State::operator=(x);
        SetType(JointOrbitClockState::TYPE);
      }
      return *this;
    }

    /// @brief Size of the orbit (Cartesian) portion of the state.
    /// @return Always 6.
    int GetOrbitStateSize() const { return ORBIT_STATE_SIZE; }

    /// @brief Size of the clock portion of the state.
    /// @return 2 or 3.
    int GetClockStateSize() const { return clock_state_size_; }

    /// @brief Extract the Cartesian orbit sub-state.
    /// @return `Cart6` containing [r; v] in this state's frame.
    Cart6 GetOrbitState() const {
      Cart6 orbit(head<ORBIT_STATE_SIZE>(), GetFrame());
      return orbit;
    }

    /// @brief Extract the clock sub-state.
    /// @return `State` of type "ClockState" containing [b, d, (dr)] with
    ///         the corresponding units, in this state's frame.
    State GetClockState() const {
      State clock(clock_state_size_);
      clock.head(clock_state_size_) = tail(clock_state_size_);
      clock.SetType("ClockState");
      clock.SetNames(clock_state_size_ == 2 ? std::vector<std::string>{"b", "d"}
                                            : std::vector<std::string>{"b", "d", "dr"});
      std::vector<std::string> joint_units = GetUnits();
      std::vector<std::string> units(joint_units.begin() + ORBIT_STATE_SIZE, joint_units.end());
      clock.SetUnits(units);
      clock.SetFrame(GetFrame());
      return clock;
    }

    /// @brief Overwrite the orbit sub-state (and this state's frame) from
    /// `orbit_state`, leaving the clock sub-state unchanged.
    /// @param orbit_state  6-element Cartesian [r; v] state.
    void SetOrbitState(const State& orbit_state) {
      LUPNT_CHECK(orbit_state.size() == ORBIT_STATE_SIZE, "Orbit state size must be 6",
                  "JointOrbitClockState");
      head(ORBIT_STATE_SIZE) = orbit_state.head(ORBIT_STATE_SIZE);
      SetFrame(orbit_state.GetFrame());
    }

    /// @brief Overwrite the clock sub-state (and its units) from
    /// `clock_state`, leaving the orbit sub-state unchanged.
    /// @param clock_state  Clock state matching this state's clock size
    ///                      (`GetClockStateSize()`).
    void SetClockState(const State& clock_state) {
      LUPNT_CHECK(clock_state.size() == clock_state_size_, "Clock state size mismatch",
                  "JointOrbitClockState");
      tail(clock_state_size_) = clock_state.head(clock_state_size_);
      std::vector<std::string> units = GetUnits();
      std::vector<std::string> clock_units = GetClockUnitsFromState(clock_state);
      std::copy(clock_units.begin(), clock_units.end(), units.begin() + ORBIT_STATE_SIZE);
      SetUnits(units);
    }

    /// @brief Mutable view of the position sub-vector. @return Position [m].
    Eigen::Ref<Vec3> r() { return head(3); }
    /// @brief Mutable view of the velocity sub-vector. @return Velocity [m/s].
    Eigen::Ref<Vec3> v() { return segment(3, 3); }
    /// @brief Mutable reference to the clock bias. @return Clock bias [s].
    Real& b() { return this->coeffRef(ORBIT_STATE_SIZE); }
    /// @brief Mutable reference to the clock drift. @return Clock drift [s/s].
    Real& d() { return this->coeffRef(ORBIT_STATE_SIZE + 1); }

    /// @brief Mutable reference to the clock drift-rate. Only valid when
    /// `GetClockStateSize() == 3`.
    /// @return Clock drift-rate [s/s^2].
    Real& dr() {
      LUPNT_CHECK(clock_state_size_ == 3, "Clock state does not contain drift rate",
                  "JointOrbitClockState");
      return this->coeffRef(ORBIT_STATE_SIZE + 2);
    }

    /// @brief Position sub-vector. @return Position [m].
    Vec3 r() const { return head(3); }
    /// @brief Velocity sub-vector. @return Velocity [m/s].
    Vec3 v() const { return segment(3, 3); }
    /// @brief Clock bias. @return Clock bias [s].
    Real b() const { return this->coeff(ORBIT_STATE_SIZE); }
    /// @brief Clock drift. @return Clock drift [s/s].
    Real d() const { return this->coeff(ORBIT_STATE_SIZE + 1); }

    /// @brief Clock drift-rate. Only valid when `GetClockStateSize() == 3`.
    /// @return Clock drift-rate [s/s^2].
    Real dr() const {
      LUPNT_CHECK(clock_state_size_ == 3, "Clock state does not contain drift rate",
                  "JointOrbitClockState");
      return this->coeff(ORBIT_STATE_SIZE + 2);
    }
  };

}  // namespace lupnt

template <> struct fmt::formatter<lupnt::State> : fmt::formatter<std::string> {
  template <typename FormatContext>
  auto format(const lupnt::State& state, FormatContext& ctx) const {
    std::ostringstream oss;
    oss << state.GetType() << "(" << state.transpose().format(lupnt::FMT_COMPACT) << ", "
        << enum_name(state.GetFrame()) << ")";
    return fmt::formatter<std::string>::format(oss.str(), ctx);
  }
};

template <> struct fmt::formatter<lupnt::Cart6> : fmt::formatter<std::string> {
  template <typename FormatContext>
  auto format(const lupnt::Cart6& state, FormatContext& ctx) const {
    std::ostringstream oss;
    oss << state.GetType() << "(" << state.transpose().format(lupnt::FMT_COMPACT) << ", "
        << enum_name(state.GetFrame()) << ")";
    return formatter<std::string>::format(oss.str(), ctx);
  }
};
