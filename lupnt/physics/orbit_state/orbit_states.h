/**
 * @file orbit_state.h
 * @author Stanford NAV LAB
 * @brief  Orbit states
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <stdexcept>

#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/state.h"

#define kOrbitStateSize 6

#define GETSET_ELEM(name, idx)                       \
  inline Real name() const { return GetVec()(idx); } \
  inline void Set_##name(Real val) { SetValue(val, idx); }

namespace lupnt {

enum class OrbitStateRepres {
  CARTESIAN = 0,
  CLASSICAL_OE,
  QUASI_NONSINGULAR_OE,
  SINGULAR_ROE,
  NONSINGULAR_OE,
  EQUINOCTIAL_OE,
  DELAUNAY_OE,
  ABSOLUTE_RELATIVE_SEPARATOR,
  RTN,
  QUASINONSINGULAR_ROE,
};

// Base class for orbit states
class OrbitState : public IState {
 private:
  Vec6 x_;
  OrbitStateRepres repres_;
  Frame frame_;
  std::vector<std::string> names_;
  std::vector<std::string> units_;

 public:
  OrbitState(const Vec6 &x, Frame coord, OrbitStateRepres repres,
             const std::vector<std::string> &names,
             const std::vector<std::string> &units)
      : x_(x), frame_(coord), repres_(repres), names_(names), units_(units) {}

  Vec6 GetVec() const { return x_; }

  VecX GetVecX() const {
    VecX x(kOrbitStateSize);
    for (int i = 0; i < kOrbitStateSize; i++) {
      x(i) = x_(i);
    }
    return x;
  }

  inline std::vector<std::string> GetNames() const { return names_; }
  inline std::vector<std::string> GetUnits() const { return units_; }
  inline int GetSize() const { return kOrbitStateSize; }
  inline Real GetValue(int i) const { return x_(i); }
  inline OrbitStateRepres GetOrbitStateRepres() const { return repres_; }
  inline Frame GetCoordSystem() const { return frame_; }

  inline void SetValue(Real val, int idx) { x_(idx) = val; }
  inline void SetVec(const Vec6 &x) { x_ = x; }

  inline void SetVecX(const VecX &x) {
    if (x.size() != kOrbitStateSize) {
      throw std::invalid_argument("Vec size does not match");
    }
    for (int i = 0; i < kOrbitStateSize; i++) {
      x_(i) = x(i);
    }
  }

  inline void SetOrbitStateRepres(const OrbitStateRepres rep) { repres_ = rep; }
  inline void SetCoordSystem(Frame sys) { frame_ = sys; }

  inline Real operator()(int idx) const { return x_(idx); }
};

// Extends OrbitState to represent an orbit in Cartesian coordinates
// The state vector is
// - $r_x$ [km]
// - $r_y$ [km]
// - $r_z$ [km]
// - $v_x$ [km/s]
// - $v_y$ [km/s]
// - $v_z$ [km/s]
class CartesianOrbitState : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {"rx", "ry", "rz",
                                                         "vx", "vy", "vz"};
  inline static const std::vector<std::string> units_ = {
      "km", "km", "km", "km/s", "km/s", "km/s"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::CARTESIAN;

 public:
  CartesianOrbitState(const Vec6 &x, Frame sys = Frame::MOON_CI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  inline Vec3 r() const { return GetVec().head(3); }
  inline Vec3 v() const { return GetVec().tail(3); }
  inline void Set_r(const Vec3 &r) {
    for (int i = 0; i < 3; i++) SetValue(r(i), i);
  }
  inline void Set_v(const Vec3 &v) {
    for (int i = 0; i < 3; i++) SetValue(v(i), i + 3);
  }
};

// Extends OrbitState to represent an orbit in classical orbital elements
// The state vector is
// - $a$ [km] (semi-major axis)
// - $e$ [-] (eccentricity)
// - $i$ [rad] (inclination)
// - $\Omega$ [rad] (right ascension of the ascending node)
// - $\omega$ [rad] (argument of periapsis)
// - $M$ [rad] (mean anomaly)
// Singular at $e \in {0,1}$, and $i \in {0, \pi}$
class ClassicalOE : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {"a",     "e", "i",
                                                         "Omega", "w", "M"};
  inline static const std::vector<std::string> units_ = {"km",  "-",   "rad",
                                                         "rad", "rad", "rad"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::CLASSICAL_OE;
  static Vec6 to_deg(const Vec6 &x, bool deg) {
    if (!deg) return x;
    Vec6 x_deg = x;
    x_deg.segment(2, 4) *= RAD;
    return x_deg;
  }

 public:
  ClassicalOE(const Vec6 &x, const Frame sys = Frame::MOON_CI, bool deg = false)
      : OrbitState(to_deg(x, deg), sys, repres_, names_, units_) {}

  GETSET_ELEM(a, 0);
  GETSET_ELEM(e, 1);
  GETSET_ELEM(i, 2);
  GETSET_ELEM(Omega, 3);
  GETSET_ELEM(w, 4);
  GETSET_ELEM(M, 5);
};

// Extends OrbitState to represent an orbit in quasi-nonsingular orbital
// elements.
// The state vector is
// - $a$ [km] (semi-major axis)
// - $u$ [-] (mean argument of latitude)
// - $e_x$ [-] (eccentricity x-component)
// - $e_y$ [-] (eccentricity y-component)
// - $i$ [rad] (inclination)
// - $\Omega$ [rad] (right ascension of the ascending node)
class QuasiNonsingOE : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {"a",  "u", "ex",
                                                         "ey", "i", "Omega"};
  inline static const std::vector<std::string> units_ = {"km", "-",   "-",
                                                         "-",  "rad", "rad"};
  static constexpr OrbitStateRepres repres_ =
      OrbitStateRepres::QUASI_NONSINGULAR_OE;

 public:
  QuasiNonsingOE(const Vec6 &x, const Frame sys = Frame::MOON_CI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(a, 0);
  GETSET_ELEM(u, 1);
  GETSET_ELEM(ex, 2);
  GETSET_ELEM(ey, 3);
  GETSET_ELEM(i, 4);
  GETSET_ELEM(Omega, 5);
};

// Extends OrbitState to represent an orbit in Delaunay orbital elements
// The state vector is
// - $l$ [rad] (mean longitude)
// - $g$ [rad] (longitude of periapsis)
// - $h$ [rad] (longitude of ascending node)
// - $L$ [rad]
// - $G$ [rad]
// - $H$ [rad]
class DelaunayOE : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {"l", "g", "h",
                                                         "L", "G", "H"};
  inline static const std::vector<std::string> units_ = {"rad", "rad", "rad",
                                                         "rad", "rad", "rad"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::DELAUNAY_OE;

 public:
  DelaunayOE(const Vec6 &x, const Frame sys = Frame::MOON_CI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(l, 0);
  GETSET_ELEM(g, 1);
  GETSET_ELEM(h, 2);
  GETSET_ELEM(L, 3);
  GETSET_ELEM(G, 4);
  GETSET_ELEM(H, 5);
};

// Extends OrbitState to represent an orbit in equinoctial orbital elements
// The state vector is
// - $a$ [km] (semi-major axis)
// - $h$ [-]
// - $k$ [-]
// - $p$ [-]
// - $q$ [-]
// - $\lambda$ [rad]
class EquinoctialOE : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {"a", "h", "k",
                                                         "p", "q", "lambda"};
  inline static const std::vector<std::string> units_ = {"km", "-", "-",
                                                         "-",  "-", "rad"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::EQUINOCTIAL_OE;

 public:
  EquinoctialOE(const Vec6 &x, const Frame sys = Frame::MOON_CI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(a, 0);
  GETSET_ELEM(h, 1);
  GETSET_ELEM(k, 2);
  GETSET_ELEM(p, 3);
  GETSET_ELEM(q, 4);
  GETSET_ELEM(lon, 5);
};

// Extends OrbitState to represent an orbit in singular relative orbital
// elements The state vector is
// - $a\delta a$ [m] (semi-major axis)
// - $a\delta M$ [m] (mean anomaly)
// - $a\delta e$ [m] (eccentricity)
// - $a\delta \omega$ [m] (argument of periapsis)
// - $a\delta i$ [m] (inclination)
// - $a\delta \Omega$ [m] (right ascension of the ascending node)
class SingularROE : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {
      "ada", "adM", "ade", "adw", "adi", "adOmega"};
  inline static const std::vector<std::string> units_ = {"m", "m", "m",
                                                         "m", "m", "m"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::SINGULAR_ROE;

 public:
  SingularROE(const Vec6 &x, const Frame sys = Frame::MOON_CI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(ada, 0);
  GETSET_ELEM(adM, 1);
  GETSET_ELEM(ade, 2);
  GETSET_ELEM(adw, 3);
  GETSET_ELEM(adi, 4);
  GETSET_ELEM(adOmega, 5);
};

// Extends OrbitState to represent an orbit in quasi-nonsingular relative
// orbital elements The state vector is
// - $a\delta a$ [m] (semi-major axis)
// - $a\delta l$ [m] (mean longitude)
// - $a\delta e_x$ [m] (eccentricity x-component)
// - $a\delta e_y$ [m] (eccentricity y-component)
// - $a\delta i_x$ [m] (inclination x-component)
// - $a\delta i_y$ [m] (inclination y-component)
class QuasiNonsingROE : public OrbitState {
 private:
  inline static const std::vector<std::string> names_ = {
      "ada", "adl", "adex", "adey", "adix", "adiy"};
  inline static const std::vector<std::string> units_ = {"m", "m", "m",
                                                         "m", "m", "m"};
  static constexpr OrbitStateRepres repres_ =
      OrbitStateRepres::QUASINONSINGULAR_ROE;

 public:
  // ada, adl, adex, adey, adix, adiy
  QuasiNonsingROE(const Vec6 &x, const Frame sys = Frame::MOON_CI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(ada, 0);
  GETSET_ELEM(adl, 1);
  GETSET_ELEM(adex, 2);
  GETSET_ELEM(adey, 3);
  GETSET_ELEM(adix, 4);
  GETSET_ELEM(adiy, 5);
};

}  // namespace lupnt