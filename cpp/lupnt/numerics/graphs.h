#pragma once

#include <functional>
#include <map>
#include <queue>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"

namespace lupnt {
  /// @brief Find the shortest sequence of representation/frame nodes connecting `start` to
  /// `end` via a breadth-first search over the given adjacency map.
  ///
  /// Used by the generic state/frame conversion machinery (e.g.
  /// `StateConverter::Convert` in conversions/state_converter.cc) to chain together a
  /// sequence of pairwise conversion functions when no direct conversion between two
  /// state/coordinate types is registered, by searching the graph of registered
  /// pairwise conversions for a path from the source type to the target type.
  ///
  /// @param start  Starting node (e.g. source state/frame type)
  /// @param end    Target node (e.g. destination state/frame type)
  /// @param map    Adjacency map: keys are (from, to) node pairs for which a direct
  ///               conversion function exists; values are unused by the search itself
  /// @return       Ordered list of nodes from `start` to `end` (inclusive) forming the
  ///               shortest path; throws (via LUPNT_CHECK) if no path exists
  template <typename T, typename U>
  std::vector<T> FindShortestPath(const T start, const T end,
                                  const std::map<std::pair<T, T>, std::function<U>>& map) {
    std::queue<T> queue;
    std::map<T, T> predecessors;
    std::map<T, bool> visited;

    queue.push(start);
    visited[start] = true;
    predecessors[start] = start;  // Start node is its own predecessor

    while (!queue.empty()) {
      T current = queue.front();
      queue.pop();

      if (current == end) {
        // Path found, reconstruct it
        std::vector<T> path;
        for (T at = end; at != start; at = predecessors[at]) {
          path.push_back(at);
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
        return path;
      }

      // Explore neighbors
      for (const auto& entry : map) {
        const auto& [repres_from, repres_to] = entry.first;
        T neighbor = repres_to;
        if (repres_from == current && !visited[neighbor]) {
          queue.push(neighbor);
          visited[neighbor] = true;
          predecessors[neighbor] = current;
        }
      }
    }
    LUPNT_CHECK(false, "Path not found from start to end representation.", "FindShortestPath");
  }

}  // namespace lupnt
