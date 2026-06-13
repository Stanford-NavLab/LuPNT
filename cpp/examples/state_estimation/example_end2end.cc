// Skeleton for an end-to-end navigation simulation. The commented sections show
// the intended order: load YAML, set up truth/filter dynamics, schedule apps,
// and run the simulation.
#include <lupnt/lupnt.h>

int main() {
  using namespace lupnt;
  spdlog::set_level(spdlog::level::debug);

  // // Path
  // auto input_path = GetInputPath("example_end2end.yaml");
  // auto output_path = SetOutputDir("example_end2end");
  // auto cfg = GetYamlConfig(input_path);

  // // Config
  // //
  // **************************************************************************
  // double t0_utc = GregorianToTime(cfg["time"]["t0_utc"].as<std::string>());
  // double t0_tai = ConvertTime(t0_utc, Time::UTC, Time::TAI);
  // double dt = cfg["time"]["dt"].as<double>();
  // double Dt = cfg["time"]["Dt"].as<double>();
  // double tf = cfg["time"]["tf"].as<double>();

  // // Setup
  // //
  // **************************************************************************
  // auto dyn_sc_true = MakePtr<NBodyDynamics>();

  // auto dyn_sc_filter = MakePtr<NBodyDynamicsT<Real>>();

  // auto sat = MakePtr<Spacecraft>();
  // auto rover = MakePtr<Rover>();

  // // Run
  // //
  // **************************************************************************
  // Simulation::Schedule(app);
  // Simulation::Run(tf);
}
