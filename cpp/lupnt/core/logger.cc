#include "lupnt/core/logger.h"

#include <sstream>

#include "lupnt/core/data_logger.h"

namespace lupnt {

  // Static member initialization
  std::chrono::steady_clock::time_point Logger::start_time_ = std::chrono::steady_clock::now();
  Logger::LogLevel Logger::log_level_ = Logger::INFO;

  std::string Logger::FormatTime(double seconds, bool short_format) {
    int hours = static_cast<int>(seconds / 3600);
    int minutes = static_cast<int>((static_cast<int>(seconds) % 3600) / 60);
    double secs = seconds - hours * 3600 - minutes * 60;

    if (!short_format) return fmt::format("{:02}:{:02}:{:05.2f}", hours, minutes, secs);
    if (hours > 0) return fmt::format("{:02}:{:02}:{:05.2f}", hours, minutes, secs);
    if (minutes > 0) return fmt::format("{:02}:{:05.2f}", minutes, secs);
    return fmt::format("{:05.2f}", secs);
  }

  std::string Logger::GetPrefix(const std::string& name, double time_sim) {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration<double>(now - start_time_).count();
    std::string elapsed_str = FormatTime(elapsed, true);
    std::stringstream ss;
    ss << "[" << Colors::CYAN << elapsed_str << Colors::RESET << "]";
    if (time_sim >= 0.0) {
      std::string sim_str = FormatTime(time_sim, true);
      ss << "[" << Colors::CYAN << sim_str << Colors::RESET << "]";
    }
    if (!name.empty()) {
      ss << "[" << Colors::GREEN << name << Colors::RESET << "]";
    }
    return ss.str();
  }

  void Logger::Log(const std::string& message, const std::string& name, double time_sim,
                   const std::string& color, LogLevel level) {
    // Prefix
    std::string prefix = GetPrefix(name, time_sim);

    // Console
    if (level >= log_level_) {
      auto progress_bars = ProgressBar::GetProgressBars();
      int n_bars = static_cast<int>(progress_bars.size());
      if (n_bars > 0) {
        std::cout << "\033[" << n_bars << "A" << std::flush;  // Move up
        std::cout << "\r\033[2K" << std::flush;               // Clear line
      }
      std::cout << prefix << " " << color << message << Colors::RESET << std::endl;
      if (n_bars > 0) std::cout << std::string(n_bars, '\n') << std::flush;
      if (level != LogLevel::ERROR) {
        for (auto bar : progress_bars) bar->Redraw();
      }
    }

    // Remove colors from message
    auto colors = {Colors::CYAN, Colors::GREEN, Colors::YELLOW, Colors::RED,
                   Colors::BLUE, Colors::WHITE, Colors::RESET};
    std::string clean_message = prefix + " " + message;
    for (auto color_code : colors) {
      // Prevent infinite loop if an empty string (like Colors::RESET on some systems) is searched.
      if (color_code.empty()) {
        continue;
      }

      size_t pos = 0;
      while ((pos = clean_message.find(color_code, pos)) != std::string::npos) {
        clean_message.erase(pos, color_code.length());
        // Do not advance pos, next iteration correctly resumes search at the current index.
      }
    }
  }

  // Static methods
  void Logger::Debug(const std::string& message, const std::string& name, double time_sim) {
    Log(message, name, time_sim, Colors::BLUE, LogLevel::DEBUG);
  }

  void Logger::Info(const std::string& message, const std::string& name, double time_sim) {
    Log(message, name, time_sim, Colors::WHITE, LogLevel::INFO);
  }

  void Logger::Warn(const std::string& message, const std::string& name, double time_sim) {
    Log(message, name, time_sim, Colors::YELLOW, LogLevel::WARNING);
  }

  void Logger::Error(const std::string& message, const std::string& name, double time_sim) {
    Log(message, name, time_sim, Colors::RED, LogLevel::ERROR);
  }

  Ptr<ProgressBar> Logger::GetProgressBar(int total, const std::string& message,
                                          const std::string& name, double time_sim) {
    Ptr<ProgressBar> bar;
    bar = MakePtr<ProgressBar>(total);

    std::stringstream ss;
    ss << GetPrefix(name, time_sim);
    if (!message.empty()) ss << " " << message;
    bar->SetDescription(ss.str());
    return bar;
  }  // namespace lupnt

  void Logger::PlaySound() { std::cout << '\a' << std::flush; }

}  // namespace lupnt
