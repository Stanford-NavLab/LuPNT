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

  Ptr<OrbitState> ConvertOrbitStateRepresentation(const Ptr<OrbitState> &state_in,
                                                  OrbitStateRepres repres_out, Real GM);

  CartesianOrbitState ConvertOrbitStateFrame(const CartesianOrbitState state, const Real epoch,
                                             const Frame frame_out);

}  // namespace lupnt
