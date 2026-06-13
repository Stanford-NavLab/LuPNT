#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "lupnt/states/state.h"

namespace lupnt {

  extern std::map<std::pair<StateType, StateType>, std::function<State(const State &, Real)>>
      absolute_conversions;

  extern std::map<std::pair<StateType, StateType>,
                  std::function<State(const State &, const State &)>>
      relative_conversions;

  State ConvertState(const State &x, StateType type_out);
  State ConvertState(const State &x, StateType type_out, Real GM);
  State ConvertState(const State &x, const State &y, StateType type_out, Real GM);

  State ConvertState(Real t_tdb, const State &x, StateType type_out, Frame frame_out);

}  // namespace lupnt
