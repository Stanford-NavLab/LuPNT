#include "converter.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <map>
#include <queue>
#include <vector>

#include "conversions.h"
#include "orbit_states.h"

#define ABSOLUTE_CONVERSION(from, to, func)        \
  {{OrbitStateRepres::from, OrbitStateRepres::to}, \
   [](const Vec6 &x, Real mu) -> Vec6 { return func(x, mu); }}

#define RELATIVE_CONVERSION(from, to, func)        \
  {{OrbitStateRepres::from, OrbitStateRepres::to}, \
   [](const Vec6 &x, const Vec6 &y) -> Vec6 { return func(x, y); }}

namespace lupnt {

std::vector<OrbitStateRepres> FindShortestPath(
    OrbitStateRepres start, OrbitStateRepres end,
    const std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
                   std::function<Vec6(const Vec6 &, Real)>>
        &absolute_conversions) {
  std::queue<OrbitStateRepres> queue;
  std::map<OrbitStateRepres, OrbitStateRepres> predecessors;
  std::map<OrbitStateRepres, bool> visited;

  queue.push(start);
  visited[start] = true;
  predecessors[start] = start;  // Start node is its own predecessor

  while (!queue.empty()) {
    OrbitStateRepres current = queue.front();
    queue.pop();

    if (current == end) {
      // Path found, reconstruct it
      std::vector<OrbitStateRepres> path;
      for (OrbitStateRepres at = end; at != start; at = predecessors[at]) {
        path.push_back(at);
      }
      path.push_back(start);
      std::reverse(path.begin(), path.end());
      return path;
    }

    // Explore neighbors
    for (const auto &entry : absolute_conversions) {
      const auto &[repres_from, repres_to] = entry.first;
      OrbitStateRepres neighbor = repres_to;
      if (repres_from == current && !visited[neighbor]) {
        queue.push(neighbor);
        visited[neighbor] = true;
        predecessors[neighbor] = current;
      }
    }
  }

  assert(false && "Path not found from start to end representation.");
}

std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
         std::function<Vec6(const Vec6 &, Real)>>
    absolute_conversions = {
        ABSOLUTE_CONVERSION(CARTESIAN, CLASSICAL_OE, Cart2Classical),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, CARTESIAN, Classical2Cart),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, QUASI_NONSINGULAR_OE,
                            Classical2QuasiNonsing),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, EQUINOCTIAL_OE,
                            Classical2Equinoctial),
        ABSOLUTE_CONVERSION(CLASSICAL_OE, DELAUNAY_OE, Classical2Delaunay),
        ABSOLUTE_CONVERSION(QUASI_NONSINGULAR_OE, CLASSICAL_OE,
                            QuasiNonsing2Classical),
        ABSOLUTE_CONVERSION(EQUINOCTIAL_OE, CLASSICAL_OE,
                            Equinoctial2Classical),
        ABSOLUTE_CONVERSION(DELAUNAY_OE, CLASSICAL_OE, Delaunay2Classical),
};

std::map<std::pair<OrbitStateRepres, OrbitStateRepres>,
         std::function<Vec6(const Vec6 &, const Vec6 &)>>
    relative_conversions = {
        RELATIVE_CONVERSION(CARTESIAN, RTN, Inertial2Rtn),
        RELATIVE_CONVERSION(QUASINONSINGULAR_ROE, CLASSICAL_OE,
                            RelQuasiNonsing2Classical),
};

Vec6 ConvertOrbitState(const Vec6 &state_in, OrbitStateRepres repres_in,
                       OrbitStateRepres repres_out, Real mu) {
  if (repres_in == repres_out) {
    return state_in;
  }

  std::vector<OrbitStateRepres> path =
      FindShortestPath(repres_in, repres_out, absolute_conversions);

  Vec6 state = state_in;
  for (size_t i = 0; i < path.size() - 1; i++) {
    state = absolute_conversions[{path[i], path[i + 1]}](state, mu);
  }

  return state;
}

Vec6 ConvertOrbitState(const Vec6 &state_in_c, const Vec6 &state_in_d,
                       OrbitStateRepres repres_in_c,
                       OrbitStateRepres repres_in_d,
                       OrbitStateRepres repres_out, Real mu) {
  // Check case
  // - (absolute_c, absolute_d) to relative_d
  // - (absolute_c, relative_d) to absolute_d
  bool to_relative =
      repres_in_c < OrbitStateRepres::ABSOLUTE_RELATIVE_SEPARATOR &&
      repres_in_d < OrbitStateRepres::ABSOLUTE_RELATIVE_SEPARATOR;

  for (const auto &entry : relative_conversions) {
    const auto &[repres_from, repres_to] = entry.first;
    auto func = entry.second;

    if (to_relative && repres_to == repres_out) {
      // state_in_c is absolute
      // state_in_d is absolute
      auto state_abs_c =
          ConvertOrbitState(state_in_c, repres_in_c, repres_from, mu);
      auto state_abs_d =
          ConvertOrbitState(state_in_d, repres_in_d, repres_from, mu);
      return func(state_abs_c, state_abs_d);
    } else if (!to_relative && repres_from == repres_in_d) {
      // state_in_c is absolute
      // state_in_d is relative
      auto state_abs_c =
          ConvertOrbitState(state_in_c, repres_in_c, repres_to, mu);
      auto state_abs_d = func(state_abs_c, state_in_d);
      return ConvertOrbitState(state_abs_d, repres_to, repres_out, mu);
    }
  }

  assert(false && "Relative conversion not found for the given input.");
}

std::shared_ptr<OrbitState> ConvertOrbitStateRepresentation(
    const std::shared_ptr<OrbitState> &state_in, OrbitStateRepres repres_out,
    Real mu) {
  Vec6 state_out = ConvertOrbitState(
      state_in->GetVec(), state_in->GetOrbitStateRepres(), repres_out, mu);
  return std::make_shared<OrbitState>(state_out, state_in->GetCoordSystem(),
                                      repres_out, state_in->GetNames(),
                                      state_in->GetUnits());
}

}  // namespace lupnt