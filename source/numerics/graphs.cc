#include "lupnt/numerics/graphs.h"

#include <algorithm>
#include <cassert>
#include <queue>
#include <vector>

namespace lupnt {

  template <typename T, typename U>
  std::vector<T> FindShortestPath(const T& start, const T& end,
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

    assert(false && "Path not found from start to end representation.");
    return {};
  }

  template std::vector<std::string> FindShortestPath<std::string, Real(Real)>(
      const std::string& start, const std::string& end,
      const std::map<std::pair<std::string, std::string>, std::function<Real(Real)>>& map);

  enum class Frame;
  template std::vector<Frame> FindShortestPath<Frame, Vec6(Real, const Vec6&)>(
      const Frame& start, const Frame& end,
      const std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6&)>>& map);

}  // namespace lupnt