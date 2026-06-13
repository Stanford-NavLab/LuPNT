#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "lupnt/core/definitions.h"

namespace lupnt {

  class ProgressBar {
  public:
    ProgressBar(int total);

    void SetWidth(int bar_width) { bar_width_ = bar_width; }
    void SetMaxUpdateFreq(double max_update_freq) { max_update_freq_ = max_update_freq; }
    void SetLeave(bool leave) { leave_ = leave; }
    void SetDescription(const std::string& description);

    void Update();
    void Update(int value);
    void Redraw();
    void Finish();
    bool IsDone() const;
    void Reset();
    ~ProgressBar();
    static std::vector<ProgressBar*> GetProgressBars() { return progress_bars_; }

  private:
    int total_;
    std::string description_;
    int row_;
    std::chrono::time_point<std::chrono::system_clock> start_time_, last_update_;
    static std::vector<ProgressBar*> progress_bars_;

    bool cursor_hidden_ = false;
    bool leave_ = false;
    int bar_width_ = 10;
    double max_update_freq_ = 0.2;
    int current_value_ = 0;
    int current_progress_ = 0;
    int value_at_last_update_ = 0;

    void Display(int value, double speed, int remaining_time, double elapsed);
    std::string FormatTime(int seconds);
  };

}  // namespace lupnt
