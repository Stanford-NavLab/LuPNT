#include <lupnt/numerics/graphs.h>

#include <algorithm>
#include <cassert>
#include <queue>
#include <vector>

namespace lupnt {

  template std::vector<std::string> FindShortestPath<std::string, Real(Real)>(
      const std::string& start, const std::string& end,
      const std::map<std::pair<std::string, std::string>, std::function<Real(Real)>>& map);

  enum class Frame;
  template std::vector<Frame> FindShortestPath<Frame, Vec6(Real, const Vec6&)>(
      const Frame& start, const Frame& end,
      const std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6&)>>& map);

}  // namespace lupnt
