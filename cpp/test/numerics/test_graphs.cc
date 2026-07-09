#include <lupnt/numerics/graphs.h>

#include <catch2/catch_test_macros.hpp>
#include <functional>
#include <map>
#include <utility>

using namespace lupnt;

namespace {

  using EdgeMap = std::map<std::pair<int, int>, std::function<void()>>;

  EdgeMap MakeEdges(const std::vector<std::pair<int, int>>& edges) {
    EdgeMap m;
    for (const auto& e : edges) m[e] = []() {};
    return m;
  }

  TEST_CASE("numerics.graphs.find_shortest_path") {
    // Chain 0 -> 1 -> 2 -> 3
    EdgeMap chain = MakeEdges({{0, 1}, {1, 2}, {2, 3}});
    auto path = FindShortestPath(0, 3, chain);
    REQUIRE(path == std::vector<int>({0, 1, 2, 3}));

    // Direct edge is preferred over a longer detour
    EdgeMap withShortcut = MakeEdges({{0, 1}, {1, 2}, {0, 2}});
    auto path2 = FindShortestPath(0, 2, withShortcut);
    REQUIRE(path2 == std::vector<int>({0, 2}));

    // Start == end yields a single-node path
    auto path3 = FindShortestPath(0, 0, chain);
    REQUIRE(path3 == std::vector<int>({0}));
  }

  TEST_CASE("numerics.graphs.no_path_throws") {
    // 3 is unreachable from 0 (edges only lead away, never to 3)
    EdgeMap disconnected = MakeEdges({{0, 1}, {1, 2}});
    REQUIRE_THROWS(FindShortestPath(0, 3, disconnected));
  }

}  // namespace
