#pragma once

#include <iostream>
#include <mutex>
#include <random>

namespace lupnt {
  class RandomEngine {
  public:
    static void SetSeed(unsigned int new_seed);
    static std::mt19937& Get();

  private:
    static std::mt19937& GetRng();
  };

}  // namespace lupnt
