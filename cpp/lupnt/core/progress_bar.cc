#include "lupnt/core/progress_bar.h"

#include "lupnt/core/error.h"

namespace lupnt {

  std::vector<ProgressBar*> ProgressBar::progress_bars_;

  ProgressBar::ProgressBar(int total)
      : total_(total), start_time_(std::chrono::system_clock::now()), last_update_(start_time_) {
    row_ = static_cast<int>(progress_bars_.size());
    std::cout << std::endl;
    Redraw();
    for (auto bar : progress_bars_) bar->Redraw();
    progress_bars_.push_back(this);
  }

  void ProgressBar::SetDescription(const std::string& description) {
    description_ = description + " ";
    Redraw();
  }

  void ProgressBar::Reset() {
    current_value_ = 0;
    current_progress_ = 0;
    value_at_last_update_ = 0;
    start_time_ = std::chrono::system_clock::now();
    last_update_ = start_time_;
    Redraw();
  }

  void ProgressBar::Update() { Update(current_value_ + 1); }

  void ProgressBar::Update(int value) {
    LUPNT_CHECK(value >= 0 && value <= total_, "Value out of range", "ProgressBar");
    cursor_hidden_ = true;
    current_value_ = value;
    auto now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_since_last_update = now - last_update_;
    int new_progress = static_cast<int>(static_cast<double>(value) / total_ * 100);

    if (new_progress > current_progress_ && elapsed_since_last_update.count() >= max_update_freq_) {
      std::chrono::duration<double> elapsed_since_start = now - start_time_;
      double iters_per_sec = value / elapsed_since_start.count();
      int remaining_time
          = iters_per_sec > 0 ? static_cast<int>((total_ - value) / iters_per_sec) : 0;

      current_progress_ = new_progress;
      last_update_ = now;
      value_at_last_update_ = value;

      Display(value, iters_per_sec, remaining_time, elapsed_since_start.count());
    }

    if (value == total_ && !leave_) Finish();
  }

  void ProgressBar::Redraw() {
    auto now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_since_start = now - start_time_;
    double iters_per_sec = current_value_ / elapsed_since_start.count();
    int remaining_time
        = iters_per_sec > 0 ? static_cast<int>((total_ - current_value_) / iters_per_sec) : 0;
    Display(current_value_, iters_per_sec, remaining_time, elapsed_since_start.count());
  }

  void ProgressBar::Finish() {
    auto now = std::chrono::system_clock::now();
    current_value_ = total_;
    current_progress_ = 100;
    std::chrono::duration<double> elapsed_since_start = now - start_time_;
    double total_iters_per_sec = total_ / elapsed_since_start.count();
    Display(total_, total_iters_per_sec, 0, elapsed_since_start.count());
#pragma omp critical
    {
      if (std::find(progress_bars_.begin(), progress_bars_.end(), this) != progress_bars_.end()) {
        progress_bars_.erase(std::remove(progress_bars_.begin(), progress_bars_.end(), this),
                             progress_bars_.end());
      }
    }
  }

  bool ProgressBar::IsDone() const { return current_value_ == total_; }

  ProgressBar::~ProgressBar() {
    // if (cursor_hidden_) std::cout << "\033[?25h" << std::flush;  // Show cursor
  }

  void ProgressBar::Display(int value, double speed, int remaining_time, double elapsed) {
#pragma omp critical
    {
      std::stringstream ss;
      std::cout << "\033[1A" << std::flush;  // Move up
      if (row_ > 0) {
        std::cout << "\033[" << row_ << "A" << std::flush;  // Move up
      }
      // Use ANSI escape code to clear the entire line
      std::cout << "\r\033[2K" << std::flush;  // Clear line

      // Draw bar
      if (!description_.empty()) ss << description_;
      ss << std::setw(3) << current_progress_ << "%|";
      int pos = static_cast<int>(static_cast<double>(current_value_) / static_cast<double>(total_)
                                 * static_cast<double>(bar_width_));
      // LUPNT_CHECK(pos >= 0 && pos <= bar_width_, "Position out of range", "ProgressBar");
      pos = std::min(std::max(pos, 0), bar_width_);
      for (int i = 0; i < pos; ++i) ss << "\u2588";
      for (int i = pos; i < bar_width_; ++i) ss << " ";
      ss << "| " << value << "/" << total_ << " [";
      ss << FormatTime(static_cast<int>(elapsed)) << "<";
      ss << FormatTime(remaining_time) << ", ";
      ss << std::fixed << std::setprecision(1) << speed << " it/s]";
      std::cout << ss.str() << std::endl;

      // Move down
      if (row_ > 0) {
        std::cout << "\033[" << row_ << "B" << std::flush;
      }
    }
  }

  std::string ProgressBar::FormatTime(int seconds) {
    std::stringstream ss;
    int hours = seconds / 3600;
    int minutes = (seconds % 3600) / 60;
    seconds %= 60;

    if (hours > 0) ss << hours << ":";
    ss << std::setw(2) << std::setfill('0') << minutes << ":";
    ss << std::setw(2) << std::setfill('0') << seconds;

    return ss.str();
  }

}  // namespace lupnt
