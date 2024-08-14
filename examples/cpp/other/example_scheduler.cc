/**
 * @file ExampleScheduler.cpp
 * @author Stanford NAV LAB
 * @brief Example of using the Scheduler class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <lupnt/core/event.h>
#include <lupnt/core/scheduler.h>

#include <functional>
#include <iostream>
#include <queue>

using namespace lupnt;

void eventFunction(double time) { std::cout << "Event executed at time: " << time << std::endl; }

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
