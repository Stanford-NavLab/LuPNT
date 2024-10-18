#include "lupnt/core/scheduler.h"

namespace lupnt {
  // Define the static member
  std::priority_queue<Event> Scheduler::events;
  Real Scheduler::time_ = 0.0;

}  // namespace lupnt
