#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitState(py::module& m) {
  py::class_<State>(
      m, "State",
      "Generic labeled state vector: values plus a type tag, per-element names/units, and a frame.")
      .def(py::init<>())
      .def(py::init<int>())
      .def(py::init<const State&>())
      .def_property(
          "type", [](const State& s) { return s.GetType(); },
          [](State& s, const StateType& v) { s.GetType() = v; },
          "State-type tag identifying the representation.")
      .def_property(
          "names", [](const State& s) { return s.GetNames(); },
          [](State& s, const std::vector<std::string>& v) { s.SetNames(v); },
          "Per-element variable names.")
      .def_property(
          "units", [](const State& s) { return s.GetUnits(); },
          [](State& s, const std::vector<std::string>& v) { s.SetUnits(v); },
          "Per-element unit strings.")
      .def_property(
          "frame", [](const State& s) { return s.GetFrame(); },
          [](State& s, const Frame& f) { s.SetFrame(f); },
          "Reference frame the state's coordinates are expressed in.");

  // ClockState
  py::class_<ClockState3>(m, "ClockState",
                          "Clock state: bias [s], drift [s/s], and drift-rate [s/s^2].")
      .def(py::init<>())
      .def(py::init<const State&>())
      .def_property(
          "type", [](const ClockState3& s) { return s.GetType(); },
          [](ClockState3& s, const StateType& v) { s.GetType() = v; },
          "State-type tag identifying the representation.")
      .def_property(
          "names", [](const ClockState3& s) { return s.GetNames(); },
          [](ClockState3& s, const std::vector<std::string>& v) { s.SetNames(v); },
          "Per-element variable names.")
      .def_property(
          "units", [](const ClockState3& s) { return s.GetUnits(); },
          [](ClockState3& s, const std::vector<std::string>& v) { s.SetUnits(v); },
          "Per-element unit strings.")
      .def_property(
          "frame", [](const ClockState3& s) { return s.GetFrame(); },
          [](ClockState3& s, const Frame& f) { s.SetFrame(f); },
          "Reference frame the state's coordinates are expressed in.");

  // AttitudeState
  py::class_<Attitude>(
      m, "AttitudeState",
      "Attitude state: scalar-first quaternion [-] plus body angular velocity [rad/s].")
      .def(py::init<>())
      .def(py::init<const State&>())
      .def_property(
          "type", [](const Attitude& s) { return s.GetType(); },
          [](Attitude& s, const StateType& v) { s.GetType() = v; },
          "State-type tag identifying the representation.")
      .def_property(
          "names", [](const Attitude& s) { return s.GetNames(); },
          [](Attitude& s, const std::vector<std::string>& v) { s.SetNames(v); },
          "Per-element variable names.")
      .def_property(
          "units", [](const Attitude& s) { return s.GetUnits(); },
          [](Attitude& s, const std::vector<std::string>& v) { s.SetUnits(v); },
          "Per-element unit strings.")
      .def_property(
          "frame", [](const Attitude& s) { return s.GetFrame(); },
          [](Attitude& s, const Frame& f) { s.SetFrame(f); },
          "Reference frame the state's coordinates are expressed in.");
}
