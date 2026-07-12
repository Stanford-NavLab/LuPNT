#include <lupnt/lupnt.h>
#include <pybind11/pybind11.h>

#include <filesystem>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitLogger(py::module& m) {
  auto logger_class
      = py::class_<Logger>(m, "Logger",
                           "Static leveled console logger with optional progress bars.")
            .def_static(
                "debug",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Debug(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0,
                "Log a message at DEBUG level, optionally tagged with a source name and sim time "
                "[s].")
            .def_static(
                "info",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Info(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0,
                "Log a message at INFO level, optionally tagged with a source name and sim time "
                "[s].")
            .def_static(
                "warn",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Warn(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0,
                "Log a message at WARNING level, optionally tagged with a source name and sim time "
                "[s].")
            .def_static(
                "error",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Error(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0,
                "Log a message at ERROR level, optionally tagged with a source name and sim time "
                "[s].")
            .def_static(
                "set_log_level",
                [](Logger::LogLevel level) {
                  if (level == Logger::LogLevel::DEBUG)
                    Logger::SetLogLevel(Logger::DEBUG);
                  else if (level == Logger::LogLevel::INFO)
                    Logger::SetLogLevel(Logger::INFO);
                  else if (level == Logger::LogLevel::WARNING)
                    Logger::SetLogLevel(Logger::WARNING);
                  else if (level == Logger::LogLevel::ERROR)
                    Logger::SetLogLevel(Logger::ERROR);
                  else
                    LUPNT_CHECK(false, "Unknown log level: " + std::string(enum_name(level)),
                                "Logger");
                },
                py::arg("level"), "Set the minimum log level that will be emitted.")
            .def_static(
                "get_progress_bar",
                [](int total, const std::string& message = "", const std::string& name = "",
                   double time_sim = -1.0) {
                  return Logger::GetProgressBar(total, message, name, time_sim);
                },
                py::arg("total"), py::arg("message") = "", py::arg("name") = "",
                py::arg("time_sim") = -1.0,
                "Create a ProgressBar with the given total number of steps.");

  logger_class.attr("LogLevel")
      = py::enum_<Logger::LogLevel>(logger_class, "LogLevel", "Logger verbosity level.")
            .value("DEBUG", Logger::DEBUG, "Debug and above.")
            .value("INFO", Logger::INFO, "Info and above.")
            .value("WARNING", Logger::WARNING, "Warning and above.")
            .value("ERROR", Logger::ERROR, "Errors only.");

  py::class_<ProgressBar, Ptr<ProgressBar>>(m, "ProgressBar", "Terminal progress bar.")
      .def("update", py::overload_cast<>(&ProgressBar::Update), "Advance the bar by one step.")
      .def("update", py::overload_cast<int>(&ProgressBar::Update), py::arg("value"),
           "Set the bar's current progress to an absolute value.")
      .def("set_description", &ProgressBar::SetDescription, py::arg("description"),
           "Set the text shown alongside the bar.")
      .def("finish", &ProgressBar::Finish, "Complete and finalize the bar.")
      .def("reset", &ProgressBar::Reset, "Reset progress back to zero.")
      .def("is_done", &ProgressBar::IsDone, "True once the bar has reached its total.")
      .def("set_width", &ProgressBar::SetWidth, py::arg("width"),
           "Set the bar width in characters.")
      .def("set_max_update_freq", &ProgressBar::SetMaxUpdateFreq, py::arg("max_update_freq"),
           "Cap the redraw frequency [Hz].")
      .def("set_leave", &ProgressBar::SetLeave, py::arg("leave"),
           "Whether to leave the finished bar on screen.");
}
