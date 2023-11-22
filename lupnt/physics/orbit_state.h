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

#include <boost/preprocessor.hpp>
#include <stdexcept>

#include "coord_converter.h"
#include "state.h"

#define GETSET_ELEM(name, idx)                          \
  inline real name() const { return GetVector()(idx); } \
  inline void Set_##name(real val) { SetValue(val, idx); }

namespace lupnt {

enum class OrbitStateRepres {
  CARTESIAN = 0,
  CLASSICAL_OE,
  QUASI_NONSINGULAR_OE,
  NONSINGULAR_OE,
  EQUINOCTIAL_OE,
  SINGULAR_ROE,
  QUASINONSINGULAR_ROE,
  DELAUNAY_OE,
};

class OrbitState;
std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
    std::shared_ptr<OrbitState> fromOrbitState, OrbitStateRepres toRepres,
    double mu = 0.0);

class OrbitState : public IState {
 private:
  Vector6real x_;
  OrbitStateRepres state_repres_;
  CoordSystem coord_sys_;

 public:
  int state_size = 6;
  OrbitState(const Vector6real &x, const CoordSystem sys,
             const OrbitStateRepres rep)
      : x_(x), coord_sys_(sys), state_repres_(rep){};
  virtual ~OrbitState(){};

  Vector6real GetVector() const { return x_; }
  inline int GetStateSize() const { return 6; }
  inline real GetValue(int i) const { return x_(i); }
  inline OrbitStateRepres GetOrbitStateRepres() const { return state_repres_; }
  inline CoordSystem GetCoordSystem() const { return coord_sys_; }

  inline void SetValue(real val, int idx) { x_(idx) = val; }
  inline void SetVector(const Vector6real &x) { x_ = x; }
  inline void SetOrbitStateRepres(const OrbitStateRepres rep) {
    state_repres_ = rep;
  }
  inline void SetCoordSystem(CoordSystem sys) { coord_sys_ = sys; }

  inline real operator()(int idx) const { return x_(idx); }

  virtual void Print(bool deg = true) const = 0;
  virtual std::shared_ptr<OrbitState> Clone() const = 0;
};

class CartesianOrbitState : public OrbitState {
 public:
  // rx, ry, rz, vx, vy, vz
  CartesianOrbitState(const Vector6real &x, CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::CARTESIAN) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<CartesianOrbitState>(*this);
  }

  inline Vector3real r() const { return GetVector().head(3); }
  inline Vector3real v() const { return GetVector().tail(3); }
  inline void Set_r(const Vector3real &r) {
    for (int i = 0; i < 3; i++) SetValue(r(i), i);
  }
  inline void Set_v(const Vector3real &v) {
    for (int i = 0; i < 3; i++) SetValue(v(i), i + 3);
  }
};

class ClassicalOE : public OrbitState {
 public:
  // a, e, i, Omega, w, M;
  ClassicalOE(const Vector6real &x, const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::CLASSICAL_OE) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<ClassicalOE>(*this);
  }

  GETSET_ELEM(a, 0);
  GETSET_ELEM(e, 1);
  GETSET_ELEM(i, 2);
  GETSET_ELEM(Omega, 3);
  GETSET_ELEM(w, 4);
  GETSET_ELEM(M, 5);
};

class QuasiNonsingularOE : public OrbitState {
 public:
  // a, u, ex, ey, i, Omega;
  QuasiNonsingularOE(const Vector6real &x,
                     const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::QUASI_NONSINGULAR_OE) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<QuasiNonsingularOE>(*this);
  }

  GETSET_ELEM(a, 0);
  GETSET_ELEM(u, 1);
  GETSET_ELEM(ex, 2);
  GETSET_ELEM(ey, 3);
  GETSET_ELEM(i, 4);
  GETSET_ELEM(Omega, 5);
};

class NonsingularOE : public OrbitState {
 public:
  // a, e1, e2, e3, e4, e5;
  NonsingularOE(const Vector6real &x, const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::NONSINGULAR_OE) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<NonsingularOE>(*this);
  }

  GETSET_ELEM(a, 0);
  GETSET_ELEM(e1, 1);
  GETSET_ELEM(e2, 2);
  GETSET_ELEM(e3, 3);
  GETSET_ELEM(e4, 4);
  GETSET_ELEM(e5, 5);
};

class DelaunayOE : public OrbitState {
 public:
  // l, g, h, L, G, H
  DelaunayOE(const Vector6real &x, const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::DELAUNAY_OE) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<DelaunayOE>(*this);
  }

  GETSET_ELEM(l, 0);
  GETSET_ELEM(g, 1);
  GETSET_ELEM(h, 2);
  GETSET_ELEM(L, 3);
  GETSET_ELEM(G, 4);
  GETSET_ELEM(H, 5);
};

class EquinoctialOE : public OrbitState {
 public:
  // a, h, k, p, q, lon
  EquinoctialOE(const Vector6real &x, const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::EQUINOCTIAL_OE) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<EquinoctialOE>(*this);
  }

  GETSET_ELEM(a, 0);
  GETSET_ELEM(h, 1);
  GETSET_ELEM(k, 2);
  GETSET_ELEM(p, 3);
  GETSET_ELEM(q, 4);
  GETSET_ELEM(lon, 5);
};

class SingularROE : public OrbitState {
 public:
  // da, dM, de, dw, di, dOmega
  SingularROE(const Vector6real &x, const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::SINGULAR_ROE) {}
  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<SingularROE>(*this);
  }

  GETSET_ELEM(da, 0);
  GETSET_ELEM(dM, 1);
  GETSET_ELEM(de, 2);
  GETSET_ELEM(dw, 3);
  GETSET_ELEM(di, 4);
  GETSET_ELEM(dOmega, 5);
};

class QuasiNonsingularROE : public OrbitState {
 public:
  // da, dl, dex, dey, dix, diy
  QuasiNonsingularROE(const Vector6real &x,
                      const CoordSystem sys = CoordSystem::NONE)
      : OrbitState(x, sys, OrbitStateRepres::QUASINONSINGULAR_ROE){};

  void Print(bool deg = true) const override;
  std::shared_ptr<OrbitState> Clone() const override {
    return std::make_shared<QuasiNonsingularROE>(*this);
  }

  GETSET_ELEM(da, 0);
  GETSET_ELEM(dl, 1);
  GETSET_ELEM(dex, 2);
  GETSET_ELEM(dey, 3);
  GETSET_ELEM(dix, 4);
  GETSET_ELEM(diy, 5);
};

}  // namespace lupnt