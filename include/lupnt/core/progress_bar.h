#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace lupnt {

  class ProgressBar {
  public:
    ProgressBar(int total, int barWidth = 50, double maxUpdateFrequency = 0.1)
        : total_(total),
          barWidth_(barWidth),
          maxUpdateFrequency_(maxUpdateFrequency),
          currentProgress_(-1),
          startTime_(std::chrono::system_clock::now()),
          lastUpdate_(startTime_),
          valueAtLastUpdate_(0),
          currentValue_(0) {}

    void Update() { Update(currentValue_ + 1); }
    void Update(int value) {
      currentValue_ = value;
      auto now = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsedSinceLastUpdate = now - lastUpdate_;
      std::chrono::duration<double> elapsedSinceStart = now - startTime_;

      int newProgress = static_cast<int>(static_cast<double>(value) / total_ * 100);
      if (newProgress > currentProgress_ && elapsedSinceLastUpdate.count() >= maxUpdateFrequency_) {
        double itersPerSecond = value / elapsedSinceStart.count();
        int remainingTime
            = itersPerSecond > 0 ? static_cast<int>((total_ - value) / itersPerSecond) : 0;

        currentProgress_ = newProgress;
        lastUpdate_ = now;
        valueAtLastUpdate_ = value;

        Display(value, itersPerSecond, remainingTime);
      }

      if (value == total_) {
        Done();
      }
    }

    void Done() {
      auto now = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsedSinceStart = now - startTime_;
      double totalItersPerSecond = total_ / elapsedSinceStart.count();
      Display(total_, totalItersPerSecond, 0);
      std::cout << std::endl;
      std::chrono::duration<double> elapsed = now - startTime_;
      std::cout << "Elapsed time: " << FormatTime(static_cast<int>(elapsed.count())) << std::endl;
    }

  private:
    int total_;
    int barWidth_;
    double maxUpdateFrequency_;
    int currentValue_;
    int currentProgress_;
    int valueAtLastUpdate_;
    std::chrono::time_point<std::chrono::system_clock> startTime_, lastUpdate_;

    void Display(int value, double speed, int remainingTime) {
      std::cout << "[";
      int pos = barWidth_ * currentProgress_ / 100;
      for (int i = 0; i < barWidth_; ++i) {
        if (i < pos)
          std::cout << "=";
        else if (i == pos)
          std::cout << ">";
        else
          std::cout << " ";
      }
      std::cout << "] " << currentProgress_ << "%, ";
      std::cout << std::fixed << std::setprecision(2) << speed << " it/s, ";
      std::cout << FormatTime(remainingTime) << " remaining\r";
      std::cout.flush();
    }

    std::string FormatTime(int seconds) {
      std::stringstream ss;
      int hours = seconds / 3600;
      int minutes = (seconds % 3600) / 60;
      seconds %= 60;

      if (hours > 0) ss << hours << "h ";
      if (hours > 0 || minutes > 0) ss << minutes << "m ";
      ss << seconds << "s";

      return ss.str();
    }
  };

}  // namespace lupnt