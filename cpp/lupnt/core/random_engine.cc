#include "lupnt/core/random_engine.h"

namespace lupnt {
  void RandomEngine::SetSeed(unsigned int new_seed) { GetRng().seed(new_seed); }

  std::mt19937& RandomEngine::GetRng() {
    thread_local static std::mt19937 rng;
    return rng;
  }

  std::mt19937& RandomEngine::Get() { return GetRng(); }
}  // namespace lupnt
