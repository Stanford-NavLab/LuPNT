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

#include "coord_converter.h"
#include "state.h"

#define kOrbitStateSize 6

#define GETSET_ELEM(name, idx)                          \
  inline real name() const { return GetVector()(idx); } \
  inline void Set_##name(real val) { SetValue(val, idx); }

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

/**
 * @class OrbitState
 * @brief Base class for orbit states
 */
class OrbitState : public IState {
 private:
  Vector6 x_;
  OrbitStateRepres repres_;
  CoordSystem coord_;
  std::array<const char *, kOrbitStateSize> names_;
  std::array<const char *, kOrbitStateSize> units_;

 public:
  OrbitState(const Vector6 &x, CoordSystem coord, OrbitStateRepres repres,
             const std::array<const char *, kOrbitStateSize> &names,
             const std::array<const char *, kOrbitStateSize> &units)
      : x_(x), coord_(coord), repres_(repres), names_(names), units_(units) {}

  Vector6 GetVector() const { return x_; }

  VectorX GetVectorX() const {
    VectorX x(kOrbitStateSize);
    for (int i = 0; i < kOrbitStateSize; i++) {
      x(i) = x_(i);
    }
    return x;
  }

  inline std::array<const char *, kOrbitStateSize> GetNames() const {
    return names_;
  }
  inline std::array<const char *, kOrbitStateSize> GetUnits() const {
    return units_;
  }
  inline int GetSize() const { return kOrbitStateSize; }
  inline real GetValue(int i) const { return x_(i); }
  inline OrbitStateRepres GetOrbitStateRepres() const { return repres_; }
  inline CoordSystem GetCoordSystem() const { return coord_; }

  inline void SetValue(real val, int idx) { x_(idx) = val; }
  inline void SetVector(const Vector6 &x) { x_ = x; }

  inline void SetVectorX(const VectorX &x) {
    if (x.size() != kOrbitStateSize) {
      throw std::invalid_argument("Vector size does not match");
    }
    for (int i = 0; i < kOrbitStateSize; i++) {
      x_(i) = x(i);
    }
  }

  inline void SetOrbitStateRepres(const OrbitStateRepres rep) { repres_ = rep; }
  inline void SetCoordSystem(CoordSystem sys) { coord_ = sys; }

  inline real operator()(int idx) const { return x_(idx); }
};

/**
 * @class CartesianOrbitState
 * @brief Extends OrbitState to represent an orbit in Cartesian coordinates
 * @details The state vector is
 * - $r_x$ [km]
 * - $r_y$ [km]
 * - $r_z$ [km]
 * - $v_x$ [km/s]
 * - $v_y$ [km/s]
 * - $v_z$ [km/s]
 */
class CartesianOrbitState : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "rx", "ry", "rz", "vx", "vy", "vz"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "km", "km", "km", "km/s", "km/s", "km/s"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::CARTESIAN;

 public:
  CartesianOrbitState(const Vector6 &x, CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  inline Vector3 r() const { return GetVector().head(3); }
  inline Vector3 v() const { return GetVector().tail(3); }
  inline void Set_r(const Vector3 &r) {
    for (int i = 0; i < 3; i++) SetValue(r(i), i);
  }
  inline void Set_v(const Vector3 &v) {
    for (int i = 0; i < 3; i++) SetValue(v(i), i + 3);
  }
};

/**
 * @class ClassicalOE
 * @brief Extends OrbitState to represent an orbit in classical orbital elements
 * @details The state vector is
 * - $a$ [km] (semi-major axis)
 * - $e$ [-] (eccentricity)
 * - $i$ [rad] (inclination)
 * - $\Omega$ [rad] (right ascension of the ascending node)
 * - $\omega$ [rad] (argument of periapsis)
 * - $M$ [rad] (mean anomaly)
 * Singular at $e \in {0,1}$, and $i \in {0, \pi}$
 */
class ClassicalOE : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "a", "e", "i", "Omega", "w", "M"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "km", "-", "rad", "rad", "rad", "rad"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::CLASSICAL_OE;

 public:
  ClassicalOE(const Vector6 &x, const CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(a, 0);
  GETSET_ELEM(e, 1);
  GETSET_ELEM(i, 2);
  GETSET_ELEM(Omega, 3);
  GETSET_ELEM(w, 4);
  GETSET_ELEM(M, 5);
};

/**
 * @class QuasiNonsingularOE
 * @brief Extends OrbitState to represent an orbit in quasi-nonsingular orbital
 * elements
 * @details The state vector is
 * - $a$ [km] (semi-major axis)
 * - $u$ [-] (mean argument of latitude)
 * - $e_x$ [-] (eccentricity x-component)
 * - $e_y$ [-] (eccentricity y-component)
 * - $i$ [rad] (inclination)
 * - $\Omega$ [rad] (right ascension of the ascending node)
 */
class QuasiNonsingularOE : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "a", "u", "ex", "ey", "i", "Omega"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "km", "-", "-", "-", "rad", "rad"};
  static constexpr OrbitStateRepres repres_ =
      OrbitStateRepres::QUASI_NONSINGULAR_OE;

 public:
  QuasiNonsingularOE(const Vector6 &x, const CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(a, 0);
  GETSET_ELEM(u, 1);
  GETSET_ELEM(ex, 2);
  GETSET_ELEM(ey, 3);
  GETSET_ELEM(i, 4);
  GETSET_ELEM(Omega, 5);
};

/**
 * @class DelaunayOE
 * @brief Extends OrbitState to represent an orbit in Delaunay orbital elements
 * @details The state vector is
 * - $l$ [rad] (mean longitude)
 * - $g$ [rad] (longitude of periapsis)
 * - $h$ [rad] (longitude of ascending node)
 * - $L$ [rad]
 * - $G$ [rad]
 * - $H$ [rad]
 */
class DelaunayOE : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "l", "g", "h", "L", "G", "H"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "rad", "rad", "rad", "rad", "rad", "rad"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::DELAUNAY_OE;

 public:
  DelaunayOE(const Vector6 &x, const CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(l, 0);
  GETSET_ELEM(g, 1);
  GETSET_ELEM(h, 2);
  GETSET_ELEM(L, 3);
  GETSET_ELEM(G, 4);
  GETSET_ELEM(H, 5);
};

/**
 * @class EquinoctialOE
 * @brief Extends OrbitState to represent an orbit in equinoctial orbital
 * elements
 * @details The state vector is
 * - $a$ [km] (semi-major axis)
 * - $h$ [-]
 * - $k$ [-]
 * - $p$ [-]
 * - $q$ [-]
 * - $\lambda$ [rad]
 */
class EquinoctialOE : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "a", "h", "k", "p", "q", "lambda"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "km", "-", "-", "-", "-", "rad"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::EQUINOCTIAL_OE;

 public:
  EquinoctialOE(const Vector6 &x, const CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(a, 0);
  GETSET_ELEM(h, 1);
  GETSET_ELEM(k, 2);
  GETSET_ELEM(p, 3);
  GETSET_ELEM(q, 4);
  GETSET_ELEM(lon, 5);
};

/**
 * @class SingularROE
 * @brief Extends OrbitState to represent an orbit in singular relative orbital
 * elements
 * @details The state vector is
 * - $a\delta a$ [m] (semi-major axis)
 * - $a\delta M$ [m] (mean anomaly)
 * - $a\delta e$ [m] (eccentricity)
 * - $a\delta \omega$ [m] (argument of periapsis)
 * - $a\delta i$ [m] (inclination)
 * - $a\delta \Omega$ [m] (right ascension of the ascending node)
 */
class SingularROE : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "ada", "adM", "ade", "adw", "adi", "adOmega"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "m", "m", "m", "m", "m", "m"};
  static constexpr OrbitStateRepres repres_ = OrbitStateRepres::SINGULAR_ROE;

 public:
  SingularROE(const Vector6 &x, const CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(ada, 0);
  GETSET_ELEM(adM, 1);
  GETSET_ELEM(ade, 2);
  GETSET_ELEM(adw, 3);
  GETSET_ELEM(adi, 4);
  GETSET_ELEM(adOmega, 5);
};

/**
 * @class QuasiNonsingularROE
 * @brief Extends OrbitState to represent an orbit in quasi-nonsingular relative
 * orbital elements
 * @details The state vector is
 * - $a\delta a$ [m] (semi-major axis)
 * - $a\delta l$ [m] (mean longitude)
 * - $a\delta e_x$ [m] (eccentricity x-component)
 * - $a\delta e_y$ [m] (eccentricity y-component)
 * - $a\delta i_x$ [m] (inclination x-component)
 * - $a\delta i_y$ [m] (inclination y-component)
 */
class QuasiNonsingularROE : public OrbitState {
 private:
  static constexpr std::array<const char *, kOrbitStateSize> names_ = {
      "ada", "adl", "adex", "adey", "adix", "adiy"};
  static constexpr std::array<const char *, kOrbitStateSize> units_ = {
      "m", "m", "m", "m", "m", "m"};
  static constexpr OrbitStateRepres repres_ =
      OrbitStateRepres::QUASINONSINGULAR_ROE;

 public:
  // ada, adl, adex, adey, adix, adiy
  QuasiNonsingularROE(const Vector6 &x, const CoordSystem sys = CoordSystem::MI)
      : OrbitState(x, sys, repres_, names_, units_) {}

  GETSET_ELEM(ada, 0);
  GETSET_ELEM(adl, 1);
  GETSET_ELEM(adex, 2);
  GETSET_ELEM(adey, 3);
  GETSET_ELEM(adix, 4);
  GETSET_ELEM(adiy, 5);
};

}  // namespace lupnt