#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "orbit_states.h"

namespace lupnt {

  extern std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
                  std::function<Vec6(const Vec6 &, Real)>>
      absolute_conversions;

  extern std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
                  std::function<Vec6(const Vec6 &, const Vec6 &)>>
      relative_conversions;

  Vec6 ConvertOrbitState(const Vec6 &state_in, OrbitStateRepres repres_in,
                         OrbitStateRepres repres_out, Real GM);

  Vec6 ConvertOrbitState(const Vec6 &state_in_c, const Vec6 &state_in_d,
                         OrbitStateRepres repres_in_c, OrbitStateRepres repres_in_d,
                         OrbitStateRepres repres_out, Real GM);

  std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
      const std::shared_ptr<OrbitState> &state_in, OrbitStateRepres repres_out, Real GM);

  static std::shared_ptr<CartesianOrbitState> ConvertOrbitStateFrame(
      const std::shared_ptr<CartesianOrbitState> state_in, const Real epoch,
      const Frame frame_out) {
    auto rv_in = state_in->GetVec();
    auto frame_in = state_in->GetCoordSystem();
    auto rv_out = ConvertFrame(epoch, rv_in, frame_in, frame_out);
    auto state_out = std::make_shared<CartesianOrbitState>(rv_out, frame_out);
    return state_out;
  }

}  // namespace lupnt
