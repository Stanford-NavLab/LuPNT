// Exercises LuPNT's logger-backed progress bars, including nested progress
// bars and interleaved log messages.
#include <thread>

#include "lupnt/core/logger.h"

using namespace lupnt;

int main() {
  auto bar1 = Logger::GetProgressBar(10, "Bar 1", "Simulation");
  for (int i = 0; i < 10; i++) {
    Logger::Info(fmt::format("i = {}", i), "Simulation", i);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    bar1->Update();
  }
  Logger::PlaySound();

  auto bar2 = Logger::GetProgressBar(10, "Bar 2", "Simulation");
  auto bar3 = Logger::GetProgressBar(10, "Bar 3", "Simulation");
  auto bar4 = Logger::GetProgressBar(10, "Bar 4", "Simulation");
  bar3->SetLeave(true);
  bar4->SetLeave(true);
  for (int i = 0; i < 10; i++) {
    Logger::Info(fmt::format("i = {}", i), "Simulation", i);
    bar3->Reset();
    for (int j = 0; j < 10; j++) {
      Logger::Info(fmt::format("j = {}", j), "Simulation", j);
      bar4->Reset();
      for (int k = 0; k < 10; k++) {
        Logger::Info(fmt::format("k = {}", k), "Simulation", k);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        bar4->Update();
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      bar3->Update();
    }
    bar2->Update();
  }
  bar2->Finish();
  bar3->Finish();
  bar4->Finish();
  Logger::PlaySound();
  return 0;
}
