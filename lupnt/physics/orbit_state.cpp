/**
 * @file OrbitState.cpp
 * @author Stanford NAV LAB
 * @brief  Orbit states
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/physics/orbit_state.h"

#include <memory>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state_utils.h"

namespace lupnt {

std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
    std::shared_ptr<OrbitState> fromOrbitState, OrbitStateRepres toRepres,
    double mu) {
  OrbitStateRepres fromRepres = fromOrbitState->GetOrbitStateRepres();
  if (fromRepres == toRepres) return fromOrbitState;
  if (fromRepres == OrbitStateRepres::CARTESIAN &&
      toRepres == OrbitStateRepres::CLASSICAL_OE) {
    if (mu == 0.0) throw std::invalid_argument("mu cannot be zero");
    std::shared_ptr<CartesianOrbitState> cs =
        std::static_pointer_cast<CartesianOrbitState>(fromOrbitState);
    return std::make_shared<ClassicalOE>(CartToCoe(*cs, mu));
  }
  if (fromRepres == OrbitStateRepres::CLASSICAL_OE &&
      toRepres == OrbitStateRepres::CARTESIAN) {
    if (mu == 0.0) throw std::invalid_argument("mu cannot be zero");
    std::shared_ptr<ClassicalOE> coe =
        std::static_pointer_cast<ClassicalOE>(fromOrbitState);
    return std::make_shared<CartesianOrbitState>(CoeToCart(*coe, mu));
  }
  return nullptr;
}

void CartesianOrbitState::Print(bool deg) const {
  std::cout << "r = " << r().segment(0, 3).transpose() << " km" << std::endl;
  std::cout << "v = " << v().segment(0, 3).transpose() << " km/s" << std::endl;
}

void ClassicalOE::Print(bool deg) const {
  if (!deg) {
    std::cout << "a     = " << a().val() << " km" << std::endl;
    std::cout << "e     = " << e().val() << std::endl;
    std::cout << "i     = " << i().val() << " rad" << std::endl;
    std::cout << "Omega = " << Omega().val() << " rad" << std::endl;
    std::cout << "w     = " << w().val() << " rad" << std::endl;
    std::cout << "M     = " << M().val() << " rad" << std::endl;
  } else {
    std::cout << "a     = " << a().val() << " km" << std::endl;
    std::cout << "e     = " << e().val() << std::endl;
    std::cout << "i     = " << i().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "Omega = " << Omega().val() * DEG_PER_RAD << " deg"
              << std::endl;
    std::cout << "w     = " << w().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "M     = " << M().val() * DEG_PER_RAD << " deg" << std::endl;
  }
}

void QuasiNonsingularOE::Print(bool deg) const {
  if (!deg) {
    std::cout << "a     = " << a().val() << " km" << std::endl;
    std::cout << "u     = " << u().val() << " rad" << std::endl;
    std::cout << "ex    = " << ex().val() << std::endl;
    std::cout << "ey    = " << ey().val() << std::endl;
    std::cout << "i     = " << i().val() << " rad" << std::endl;
    std::cout << "Omega = " << Omega().val() << " rad" << std::endl;
  } else {
    std::cout << "a     = " << a().val() << " km" << std::endl;
    std::cout << "u     = " << u().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "ex    = " << ex().val() << std::endl;
    std::cout << "ey    = " << ey().val() << std::endl;
    std::cout << "i     = " << i().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "Omega = " << Omega().val() * DEG_PER_RAD << " deg"
              << std::endl;
  }
}

void NonsingularOE::Print(bool deg) const {
  if (!deg) {
    std::cout << "a  =                     = " << a().val() << " km"
              << std::endl;
    std::cout << "e1 = E cos(Omega + w)    = " << e1().val() << " rad"
              << std::endl;
    std::cout << "e2 = E sin(Omega + w)    = " << e2().val() << " rad"
              << std::endl;
    std::cout << "e3 = sin(i/2) sin(Omega) = " << e3().val() << std::endl;
    std::cout << "e4 = sin(i/2) cos(Omega) = " << e4().val() << std::endl;
    std::cout << "e5 = Omega + w + M       = " << e5().val() << " rad"
              << std::endl;
  } else {
    std::cout << "a  =                     = " << a().val() << " km"
              << std::endl;
    std::cout << "e1 = E cos(Omega + w)    = " << e1().val() * DEG_PER_RAD
              << " deg" << std::endl;
    std::cout << "e2 = E sin(Omega + w)    = " << e2().val() * DEG_PER_RAD
              << " deg" << std::endl;
    std::cout << "e3 = sin(i/2) sin(Omega) = " << e3().val() << std::endl;
    std::cout << "e4 = sin(i/2) cos(Omega) = " << e4().val() << std::endl;
    std::cout << "e5 = Omega + w + M       = " << e5().val() * DEG_PER_RAD
              << " deg" << std::endl;
  }
}

void EquinoctialOE::Print(bool deg) const {
  if (!deg) {
    std::cout << "a   = " << a().val() << " km" << std::endl;
    std::cout << "h   = " << h().val() << std::endl;
    std::cout << "k   = " << k().val() << std::endl;
    std::cout << "p   = " << p().val() << std::endl;
    std::cout << "q   = " << q().val() << std::endl;
    std::cout << "lon = " << lon().val() << " rad" << std::endl;
  } else {
    std::cout << "a   = " << a().val() << " km" << std::endl;
    std::cout << "h   = " << h().val() << std::endl;
    std::cout << "k   = " << k().val() << std::endl;
    std::cout << "p   = " << p().val() << std::endl;
    std::cout << "q   = " << q().val() << std::endl;
    std::cout << "lon = " << lon().val() * DEG_PER_RAD << " deg" << std::endl;
  }
}

void QuasiNonsingularROE::Print(bool deg) const {
  if (!deg) {
    std::cout << "da  = " << da().val() << " km" << std::endl;
    std::cout << "dl  = " << dl().val() << " rad" << std::endl;
    std::cout << "dex = " << dex().val() << std::endl;
    std::cout << "dey = " << dey().val() << std::endl;
    std::cout << "dix = " << dix().val() << std::endl;
    std::cout << "diy = " << diy().val() << std::endl;
  } else {
    std::cout << "da  = " << da().val() << " km" << std::endl;
    std::cout << "dl  = " << dl().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "dex = " << dex().val() << std::endl;
    std::cout << "dey = " << dey().val() << std::endl;
    std::cout << "dix = " << dix().val() << std::endl;
    std::cout << "diy = " << diy().val() << std::endl;
  }
}

void DelaunayOE::Print(bool deg) const {
  if (!deg) {
    std::cout << "l = " << l().val() << " rad" << std::endl;
    std::cout << "g = " << g().val() << " rad" << std::endl;
    std::cout << "h = " << h().val() << " rad" << std::endl;
    std::cout << "L = " << L().val() << std::endl;
    std::cout << "G = " << G().val() << std::endl;
    std::cout << "H = " << H().val() << std::endl;
  } else {
    std::cout << "l = " << l().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "g = " << g().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "h = " << h().val() * DEG_PER_RAD << " deg" << std::endl;
    std::cout << "L = " << L().val() << std::endl;
    std::cout << "G = " << G().val() << std::endl;
    std::cout << "H = " << H().val() << std::endl;
  }
}

}  // namespace lupnt