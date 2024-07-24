#pragma once

#include <functional>
#include <map>
#include <queue>
#include <vector>

#include "lupnt/core/constants.h"

namespace lupnt {
template <typename T, typename U>
std::vector<T> FindShortestPath(
    const T& start, const T& end,
    const std::map<std::pair<T, T>, std::function<U>>& map);

}  // namespace lupnt