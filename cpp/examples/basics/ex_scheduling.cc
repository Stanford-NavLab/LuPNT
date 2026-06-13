// Demonstrates the event scheduler by recursively scheduling the same callback
// at a time-dependent cadence.
#include <lupnt/lupnt.h>
using namespace lupnt;

void func(Simulation& sim, Real time) {
  std::cout << time << std::endl;  // Print
  auto f = [&sim](Real next_time) { func(sim, next_time); };
  sim.Schedule(time + 1.0 / time, f);  // Reschedule
}

int main() {
  Simulation sim;
  auto f = [&sim](Real time) { func(sim, time); };
  sim.Schedule(1.2345, f);
  sim.SetDuration(3.0);
  sim.Run();
}
