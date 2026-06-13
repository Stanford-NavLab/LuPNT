#pragma once
#include <fmt/format.h>
#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <string>

#include "lupnt/core/definitions.h"
#include "lupnt/core/progress_bar.h"

static const auto LUPNT_START_TIME = std::chrono::steady_clock::now();

namespace lupnt {
  // ANSI color codes
  namespace Colors {
    const std::string RESET = "\033[0m";
    const std::string CYAN = "\033[36m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string RED = "\033[31m";
    const std::string BLUE = "\033[34m";
    const std::string WHITE = "\033[37m";
  }  // namespace Colors

  class Logger {
  public:
    enum LogLevel { DEBUG = 0, INFO = 1, WARNING = 2, ERROR = 3 };

    // Static methods for logging
    static void SetLogLevel(LogLevel level) { log_level_ = level; }
    static void Debug(const std::string& message, const std::string& name = "",
                      double time_sim = -1.0);
    static void Info(const std::string& message, const std::string& name = "",
                     double time_sim = -1.0);
    static void Warn(const std::string& message, const std::string& name = "",
                     double time_sim = -1.0);
    static void Error(const std::string& message, const std::string& name = "",
                      double time_sim = -1.0);
    static Ptr<ProgressBar> GetProgressBar(int total, const std::string& message = "",
                                           const std::string& name = "", double time_sim = -1.0);

    static void CheckProgressBars();
    static void PlaySound();

  private:
    static std::chrono::steady_clock::time_point start_time_;
    static LogLevel log_level_;
    static std::string FormatTime(double seconds, bool short_format = false);
    static std::string GetPrefix(const std::string& name = "", double time_sim = -1.0);
    static void Log(const std::string& message, const std::string& name = "",
                    double time_sim = -1.0, const std::string& color = "", LogLevel level = INFO);
  };
}  // namespace lupnt
