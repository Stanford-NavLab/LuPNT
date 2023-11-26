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

}  // namespace lupnt