#include <lupnt/core/Event.h>
#include <lupnt/core/Scheduler.h>

#include <functional>
#include <iostream>
#include <queue>

using namespace LPT;

void eventFunction(double time) {
  std::cout << "Event executed at time: " << time << std::endl;
}

int main() {
  Scheduler scheduler;

  // Schedule three events at different times
  scheduler.ScheduleEvent(Event(5.0, eventFunction));
  scheduler.ScheduleEvent(Event(2.0, eventFunction));
  scheduler.ScheduleEvent(Event(8.0, eventFunction));

  // Run the simulation
  scheduler.RunSimulation(10.0);

  return 0;
}
