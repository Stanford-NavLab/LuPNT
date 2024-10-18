#include <lupnt/lupnt.h>
using namespace lupnt;

void func(Real time, Scheduler &sch) {
  std::cout << time << ", ";                        // Print
  auto f = [&sch](Real time) { func(time, sch); };  // Function
  sch.ScheduleEvent(Event(time + 1.0 / time, f));   // Reschedule
}

int main() {
  Scheduler sch;
  auto f = [&sch](Real time) { func(time, sch); };
  sch.ScheduleEvent(Event(1.2345, f));
  sch.RunSimulation(3.0);
}
