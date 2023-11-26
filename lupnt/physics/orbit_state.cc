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
    const std::shared_ptr<OrbitState> &state_in, OrbitStateRepres repres_out,
    double mu) {
  OrbitStateRepres repres_in = state_in->GetOrbitStateRepres();

  if (repres_in == repres_out) {
    return std::make_shared<OrbitState>(*state_in);
  }

  std::shared_ptr<OrbitState> state_out = nullptr;

  if (repres_in == OrbitStateRepres::CARTESIAN &&
      repres_out == OrbitStateRepres::CLASSICAL_OE) {
    if (mu == 0.0) throw std::invalid_argument("mu cannot be zero");
    std::shared_ptr<CartesianOrbitState> cs =
        std::static_pointer_cast<CartesianOrbitState>(state_in);
    return std::make_shared<ClassicalOE>(CartToCoe(*cs, mu));
  }
  if (repres_in == OrbitStateRepres::CLASSICAL_OE &&
      repres_out == OrbitStateRepres::CARTESIAN) {
    if (mu == 0.0) throw std::invalid_argument("mu cannot be zero");
    std::shared_ptr<ClassicalOE> coe =
        std::static_pointer_cast<ClassicalOE>(state_in);
    return std::make_shared<CartesianOrbitState>(CoeToCart(*coe, mu));
  }
  return nullptr;
}

}  // namespace lupnt