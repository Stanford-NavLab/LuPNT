#include <lupnt/lupnt.h>
#include <pybind11/pybind11.h>

#include <filesystem>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitLogger(py::module& m) {
  auto logger_class
      = py::class_<Logger>(m, "Logger")
            .def_static(
                "debug",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Debug(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0)
            .def_static(
                "info",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Info(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0)
            .def_static(
                "warn",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Warn(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0)
            .def_static(
                "error",
                [](const std::string& message, const std::string& name = "",
                   double time_sim = -1.0) { Logger::Error(message, name, time_sim); },
                py::arg("message"), py::arg("name") = "", py::arg("time_sim") = -1.0)
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
                py::arg("level"))
            .def_static(
                "get_progress_bar",
                [](int total, const std::string& message = "", const std::string& name = "",
                   double time_sim = -1.0) {
                  return Logger::GetProgressBar(total, message, name, time_sim);
                },
                py::arg("total"), py::arg("message") = "", py::arg("name") = "",
                py::arg("time_sim") = -1.0);

  logger_class.attr("LogLevel") = py::enum_<Logger::LogLevel>(logger_class, "LogLevel")
                                      .value("DEBUG", Logger::DEBUG)
                                      .value("INFO", Logger::INFO)
                                      .value("WARNING", Logger::WARNING)
                                      .value("ERROR", Logger::ERROR);

  py::class_<ProgressBar, Ptr<ProgressBar>>(m, "ProgressBar")
      .def("update", py::overload_cast<>(&ProgressBar::Update))
      .def("update", py::overload_cast<int>(&ProgressBar::Update), py::arg("value"))
      .def("set_description", &ProgressBar::SetDescription, py::arg("description"))
      .def("finish", &ProgressBar::Finish)
      .def("reset", &ProgressBar::Reset)
      .def("is_done", &ProgressBar::IsDone)
      .def("set_width", &ProgressBar::SetWidth, py::arg("width"))
      .def("set_max_update_freq", &ProgressBar::SetMaxUpdateFreq, py::arg("max_update_freq"))
      .def("set_leave", &ProgressBar::SetLeave, py::arg("leave"));
}
