#include <lupnt/lupnt.h>

#include <string>
#include <vector>

#include "../utils.cc"

double ABS_TOL = 1e-6;
double REL_TOL = 1e-3;

TEST_CASE("Antenna") {
  using namespace lupnt;
  Antenna ant("Block-IIR-M_ACE");

  RequireNear(ant.ComputeGain(358 * RAD, 74 * RAD), -18.68, ABS_TOL);
  RequireNear(ant.ComputeGain(359 * RAD, 74 * RAD), -18.73, ABS_TOL);
  RequireNear(ant.ComputeGain(358.5 * RAD, 74 * RAD), 0.5 * (-18.68 - 18.73), ABS_TOL);
  RequireNear(ant.ComputeGain(359.5 * RAD, 74 * RAD), 0.5 * (-20.45 - 18.73), ABS_TOL);
  RequireNear(ant.ComputeGain(359 * RAD, 74.5 * RAD), 0.5 * (-17.48 - 18.73), ABS_TOL);

  std::vector<std::string> names
      = {"Parabora_S_d10", "Parabora_S_d100", "Block-IIA_ACE", "Block-IIR-M_ACE",
         "BEIDOU_IGSO",    "BEIDOU_MEO",      "GALLILEO",      "moongpsr",
         "DSN-S",          "DSN-X",           "LGPS",          "Patch_22_RHCP_8025MHz"};
  VecX phi = VecX::LinSpaced(-90, 90, 181) * RAD;
  VecX theta = VecX::LinSpaced(0, 360, 361) * RAD;
  for (std::string name : names) {
    Antenna ant(name);
    MatX gain_mat = ant.ComputeGain(theta, phi);
    for (int i = 0; i < gain_mat.rows(); i++) {
      for (int j = 0; j < gain_mat.cols(); j++) {
        REQUIRE(gain_mat(i, j).val() < 1e6);
        REQUIRE(gain_mat(i, j).val() > -1e6);
      }
    }
    Real gain = ant.ComputeGain(0.0, 0.0);
    REQUIRE(gain.val() < 1e6);
    REQUIRE(gain.val() > -1e6);
  }
}
