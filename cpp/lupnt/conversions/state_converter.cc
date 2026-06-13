#include "lupnt/conversions/state_converter.h"

#include <functional>
#include <map>
#include <vector>

#include "lupnt/conversions/state_conversions.h"
#include "lupnt/environment/body.h"
#include "lupnt/numerics/graphs.h"
#include "lupnt/states/state.h"

#define ABSOLUTE_CONVERSION(from, to, func) \
  {{from, to}, [](const State& x, Real GM) -> State { return func(x, GM); }}

#define RELATIVE_CONVERSION(from, to, func) \
  {{from, to}, [](const State& x, const State& y) -> State { return func(x, y); }}

namespace lupnt {

  std::map<std::pair<StateType, StateType>, std::function<State(const State&, Real)>>
      absolute_conversions = {
          ABSOLUTE_CONVERSION(Cart6::TYPE, ClassicalOE::TYPE, CartToClassical),
          ABSOLUTE_CONVERSION(ClassicalOE::TYPE, Cart6::TYPE, ClassicalToCart),
          ABSOLUTE_CONVERSION(ClassicalOE::TYPE, QuasiNonsingularOE::TYPE, ClassicalToQuasiNonsing),
          ABSOLUTE_CONVERSION(ClassicalOE::TYPE, EquinoctialOE::TYPE, ClassicalToEquinoctial),
          ABSOLUTE_CONVERSION(ClassicalOE::TYPE, DelaunayOE::TYPE, ClassicalToDelaunay),
          ABSOLUTE_CONVERSION(QuasiNonsingularOE::TYPE, ClassicalOE::TYPE, QuasiNonsingToClassical),
          ABSOLUTE_CONVERSION(EquinoctialOE::TYPE, ClassicalOE::TYPE, EquinoctialToClassical),
          ABSOLUTE_CONVERSION(DelaunayOE::TYPE, ClassicalOE::TYPE, DelaunayToClassical),
  };

  std::map<std::pair<StateType, StateType>, std::function<State(const State&, const State&)>>
      relative_conversions = {
          RELATIVE_CONVERSION(Cart6::TYPE, RelCart6::TYPE, InertialToSynodic),
          RELATIVE_CONVERSION(QuasiNonsingularOE::TYPE, ClassicalOE::TYPE,
                              RelQuasiNonsingToClassical),
  };

  State ConvertState(const State& x, StateType type_out) {
    if (x.GetType() == type_out) return x;
    return ConvertState(x, type_out, GetBodyGM(GetFrameCenter(x.GetFrame())));
  }

  State ConvertState(Real t_tdb, const State& x, StateType type_out, Frame frame_out) {
    if (x.GetType() == type_out && x.GetFrame() == frame_out) return x;
    State state = ConvertState(x, type_out);
    return ConvertFrame(t_tdb, state, frame_out);
  }

  State ConvertState(const State& x, StateType type_out, Real GM) {
    if (x.GetType() == type_out) return x;
    std::vector<StateType> path = FindShortestPath(x.GetType(), type_out, absolute_conversions);
    State state = x;
    for (size_t i = 0; i < path.size() - 1; i++)
      state = absolute_conversions[{path[i], path[i + 1]}](state, GM);
    return state;
  }

  bool IsAbsolute(const State& x) {
    StateType type = x.GetType();
    if (type == RelCart6::TYPE || SingularROE::TYPE == type || type == QuasiNonsingROE::TYPE)
      return false;
    return true;
  }

  State ConvertState(const State& x, const State& y, StateType type_out, Real GM) {
    bool y_absolute = IsAbsolute(y);

    for (const auto& entry : relative_conversions) {
      const auto& [type_1, type_2] = entry.first;
      auto func = entry.second;

      if (y_absolute && type_2 == type_out) {
        auto x_abs = ConvertState(x, type_1, GM);
        auto y_abs = ConvertState(y, type_1, GM);
        return func(x_abs, y_abs);
      } else if (!y_absolute && type_1 == y.GetType()) {
        auto x_abs = ConvertState(x, type_2, GM);
        auto y_abs = func(x_abs, y);
        return ConvertState(y_abs, type_out, GM);
      }
    }
    LUPNT_CHECK(false, "Conversion not found for the given input.", "OrbitStateConverter");
  }
}  // namespace lupnt
