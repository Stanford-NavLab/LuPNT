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

#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/state.h"

#define kOrbitStateSize 6

#define GETSET_ELEM(name, idx)                \
  Real name() const { return GetVec()(idx); } \
  void Set_##name(Real val) { SetValue(idx, val); }

namespace lupnt {

  enum class OrbitStateRepres {
    CARTESIAN,
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

  std::ostream &operator<<(std::ostream &os, const OrbitStateRepres &repres);

  // Base class for orbit states
  class OrbitState : public IState {
  private:
    Vec6 x_;
    Frame frame_;
    OrbitStateRepres repres_;
    std::array<const char *, 6> names_;
    std::array<const char *, 6> units_;

  public:
    OrbitState(const Vec6 &x, Frame coord, OrbitStateRepres repres,
               const std::array<const char *, 6> &names, const std::array<const char *, 6> &units);

    // Overrides
    int GetSize() const override;
    VecX GetVec() const override;
    void SetVec(const VecX &x) override;
    Real GetValue(int idx) const override;
    void SetValue(int idx, Real val) override;

    Vec6 GetVec6() const;
    Frame GetFrame() const;
    std::array<const char *, 6> GetNames() const;
    std::array<const char *, 6> GetUnits() const;
    StateType GetStateType() const override;
    OrbitStateRepres GetOrbitStateRepres() const;

    void SetOrbitStateRepres(const OrbitStateRepres rep);
    void SetCoordSystem(Frame frame);

    Real operator()(int idx) const;
    std::ostream &operator<<(std::ostream &os) const;
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
    static constexpr std::array<const char *, 6> names_ = {"rx", "ry", "rz", "vx", "vy", "vz"};
    static constexpr std::array<const char *, 6> units_
        = {"km", "km", "km", "km/s", "km/s", "km/s"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::CARTESIAN;

  public:
    CartesianOrbitState(const Vec6 &x, Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

    Vec3 r() const { return GetVec().head(3); }
    Vec3 v() const { return GetVec().tail(3); }
    void Set_r(const Vec3 &r) {
      for (int i = 0; i < 3; i++) SetValue(i, r(i));
    }
    void Set_v(const Vec3 &v) {
      for (int i = 0; i < 3; i++) SetValue(i + 3, v(i));
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
    static constexpr std::array<const char *, 6> names_ = {"a", "e", "i", "Omega", "w", "M"};
    static constexpr std::array<const char *, 6> units_ = {"km", "-", "rad", "rad", "rad", "rad"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::CLASSICAL_OE;

  public:
    ClassicalOE(const Vec6 &x, const Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

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
    static constexpr std::array<const char *, 6> names_ = {"a", "u", "ex", "ey", "i", "Omega"};
    static constexpr std::array<const char *, 6> units_ = {"km", "-", "-", "-", "rad", "rad"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::QUASI_NONSINGULAR_OE;

  public:
    QuasiNonsingOE(const Vec6 &x, const Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

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
    static constexpr std::array<const char *, 6> names_ = {"l", "g", "h", "L", "G", "H"};
    static constexpr std::array<const char *, 6> units_
        = {"rad", "rad", "rad", "rad", "rad", "rad"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::DELAUNAY_OE;

  public:
    DelaunayOE(const Vec6 &x, const Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

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
    static constexpr std::array<const char *, 6> names_ = {"a", "h", "k", "p", "q", "lambda"};
    static constexpr std::array<const char *, 6> units_ = {"km", "-", "-", "-", "-", "rad"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::EQUINOCTIAL_OE;

  public:
    EquinoctialOE(const Vec6 &x, const Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

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
    static constexpr std::array<const char *, 6> names_
        = {"ada", "adM", "ade", "adw", "adi", "adOmega"};
    static constexpr std::array<const char *, 6> units_ = {"m", "m", "m", "m", "m", "m"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::SINGULAR_ROE;

  public:
    SingularROE(const Vec6 &x, const Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

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
    static constexpr std::array<const char *, 6> names_
        = {"ada", "adl", "adex", "adey", "adix", "adiy"};
    static constexpr std::array<const char *, 6> units_ = {"m", "m", "m", "m", "m", "m"};
    static constexpr OrbitStateRepres repres_ = OrbitStateRepres::QUASINONSINGULAR_ROE;

  public:
    // ada, adl, adex, adey, adix, adiy
    QuasiNonsingROE(const Vec6 &x, const Frame frame = Frame::MOON_CI)
        : OrbitState(x, frame, repres_, names_, units_) {}

    GETSET_ELEM(ada, 0);
    GETSET_ELEM(adl, 1);
    GETSET_ELEM(adex, 2);
    GETSET_ELEM(adey, 3);
    GETSET_ELEM(adix, 4);
    GETSET_ELEM(adiy, 5);
  };

}  // namespace lupnt
