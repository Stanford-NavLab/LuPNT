#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace lupnt {

  class ProgressBar {
  public:
    ProgressBar(int total, int bar_width = 20, double max_update_freq = 0.2)
        : total_(total),
          bar_width_(bar_width),
          max_update_freq_(max_update_freq),
          current_value_(0),
          current_progress_(-1),
          value_at_last_update_(0),
          start_time_(std::chrono::system_clock::now()),
          last_update_(start_time_) {}

    void SetDescription(const std::string &description) { description_ = description + " "; }
    void Update() { Update(current_value_ + 1); }
    void Update(int value) {
      current_value_ = value;
      auto now = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_since_last_update = now - last_update_;
      int new_progress = static_cast<int>(static_cast<double>(value) / total_ * 100);

      if (new_progress > current_progress_
          && elapsed_since_last_update.count() >= max_update_freq_) {
        std::chrono::duration<double> elapsed_since_start = now - start_time_;
        double itersPerSecond = value / elapsed_since_start.count();
        int remaining_time
            = itersPerSecond > 0 ? static_cast<int>((total_ - value) / itersPerSecond) : 0;

        current_progress_ = new_progress;
        last_update_ = now;
        value_at_last_update_ = value;

        Display(value, itersPerSecond, remaining_time, elapsed_since_start.count());
      }

      if (value == total_) {
        Finish();
      }
    }

    void Finish() {
      auto now = std::chrono::system_clock::now();
      current_value_ = total_;
      current_progress_ = 100;
      std::chrono::duration<double> elapsed_since_start = now - start_time_;
      double total_iters_per_sec = total_ / elapsed_since_start.count();
      Display(total_, total_iters_per_sec, 0, elapsed_since_start.count());
      std::cout << std::endl;
    }

  private:
    std::string description_;
    int total_;
    int bar_width_;
    double max_update_freq_;
    int current_value_;
    int current_progress_;
    int value_at_last_update_;
    std::chrono::time_point<std::chrono::system_clock> start_time_, last_update_;

    void Display(int value, double speed, int remaining_time, double elapsed) {
      std::cout << description_;
      std::cout << "[";
      int pos = bar_width_ * current_progress_ / 100 + 1;
      for (int i = 0; i < bar_width_; ++i) {
        if (i < pos)
          std::cout << "=";
        else if (i == pos)
          std::cout << ">";
        else
          std::cout << " ";
      }
      std::cout << "] " << value << "/" << total_ << ", " << std::setw(3) << current_progress_
                << "%, ";
      std::cout << std::fixed << std::setprecision(1) << speed << " it/s, ";
      std::cout << FormatTime(static_cast<int>(elapsed)) << " elapsed, ";
      std::cout << FormatTime(remaining_time) << " left\r";
      std::cout << std::flush;
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
