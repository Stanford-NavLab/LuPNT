#pragma once

#include <numeric>
#include <vector>

namespace lupnt {

  /** @brief Helper function that returns a vector of given size and type
   *
   * @tparam T The type of element
   * @param size The size of the vector to return
   * @returns a vector of given size and type
   */
  template <typename T> std::vector<T> make_vector(std::size_t size) {
    std::vector<T> v(size, 0);
    std::iota(v.begin(), v.end(), 0);
    return v;
  }

}  // namespace lupnt
