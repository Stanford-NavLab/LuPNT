#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitAgents(py::module &m) {
  // Agent
  //   py::class_<Agent>(m, "Agent")
  //       .def(py::init<>())
  //       .def("get_epoch", &Agent::GetEpoch)
  //       .def("set_epoch", &Agent::SetEpoch, py::arg("t"))
  //       .def("get_name", &Agent::GetName)
  //       .def("set_name", &Agent::SetName, py::arg("name"))
  //       .def("get_state", &Agent::GetState)
  //       .def("set_state", &Agent::SetState, py::arg("x"))
  //       .def("get_clock", &Agent::GetClock)
  //       .def("set_clock", &Agent::SetClock, py::arg("clk"))
  //       .def("get_attitude", &Agent::GetAttitude)
  //       .def("set_attitude", &Agent::SetAttitude, py::arg("attitude"))
  //       .def("get_dynamics", &Agent::GetDynamics, py::return_value_policy::reference)
  //       .def("set_dynamics", &Agent::SetDynamics, py::arg("dyn"))
  //       .def("get_clock_dynamics", &Agent::GetClockDynamics, py::return_value_policy::reference)
  //       .def("set_clock_dynamics", &Agent::SetClockDynamics, py::arg("dyn"))
  //       .def("get_devices", &Agent::GetDevices)
  //       .def("get_transmitters", &Agent::GetTransmitters)
  //       .def("get_receivers", &Agent::GetReceivers)
  //       .def("get_transponders", &Agent::GetTransponders)
  //       .def("get_applications", &Agent::GetApplications)
  //       .def("set_body_id", &Agent::SetBodyId, py::arg("body_id"))
  //       .def("get_body_id", &Agent::GetBodyId)
  //       .def("add_application", &Agent::AddApplication, py::arg("app"))
  //       .def("add_device", &Agent::AddDevice, py::arg("device"))
  //       .def("propagate", &Agent::Propagate, py::arg("t"))
  //       .def("get_state_at", &Agent::GetStateAt, py::arg("t"))
  //       .def("get_attitude_at", &Agent::GetAttitudeAt, py::arg("t"))
  //       .def("get_clock_at", &Agent::GetClockAt, py::arg("t"))
  //       .def("log", &Agent::Log);

  //   // AgentGroup
  //   py::class_<AgentGroup>(m, "AgentGroup")
  //       .def(py::init<>())
  //       .def("set_epoch", &AgentGroup::SetEpoch, py::arg("epoch"))
  //       .def("set_channel", &AgentGroup::SetChannel, py::arg("channel"))
  //       .def("set_dynamics", &AgentGroup::SetDynamics, py::arg("dyn"))
  //       .def("get_epoch", &AgentGroup::GetEpoch)
  //       .def("get_size", &AgentGroup::GetSize)
  //       .def("get_agent", &AgentGroup::GetAgent, py::arg("i"),
  //       py::return_value_policy::reference) .def("get_channel", &AgentGroup::GetChannel,
  //       py::return_value_policy::reference) .def("get_dynamics", &AgentGroup::GetDynamics,
  //       py::return_value_policy::reference) .def("propagate", &AgentGroup::Propagate,
  //       py::arg("t"));
}
