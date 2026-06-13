// Compares LuPNT's analytic body-fixed frame orientation with the equivalent
// SPICE frame conversion for one planetary body.
#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  // Change this string to "Venus" to exercise another built-in body frame.
  std::string body = "Mars";
  BodyId id;
  Frame fixed_frame;
  std::string fixed_frame_str;

  if (body == "Mars") {
    id = BodyId::MARS;
    fixed_frame = Frame::MARS_FIXED;
    fixed_frame_str = "IAU_MARS";
  } else if (body == "Venus") {
    id = BodyId::VENUS;
    fixed_frame = Frame::VENUS_FIXED;
    fixed_frame_str = "IAU_VENUS";
  } else {
    throw std::runtime_error("Invalid body");
  }

  Real t_tdb = 0.1 * DAYS_YEAR * SECS_DAY;
  Vec4 angles = PlanetOrientation(id, t_tdb);

  std::cout << "Planet: " << body << std::endl;
  std::cout << "Epoch:  "
            << "TDB=" << t_tdb << std::endl;

  std::cout << "Orientation at t_tdb (lupnt) = " << t_tdb << " s" << std::endl;
  std::cout << "  alpha0 = " << DEG * WrapToPi(angles(0)) << std::endl;
  std::cout << "  delta0 = " << DEG * WrapToPi(angles(1)) << std::endl;
  std::cout << "  W = " << DEG * WrapToPi(angles(2)) << std::endl;
  std::cout << "  Wdot = " << DEG * angles(3) << std::endl;
  std::cout << " " << std::endl;

  Vec3d angles_spice = spice::GetPlanetOrientation(id, t_tdb);
  std::cout << "Orientation at t_tdb (SPICE) = " << t_tdb << " s" << std::endl;
  std::cout << "  alpha0 = " << DEG * WrapToPi(angles_spice(0)) << std::endl;
  std::cout << "  delta0 = " << DEG * WrapToPi(angles_spice(1)) << std::endl;
  std::cout << "  W = " << DEG * WrapToPi(angles_spice(2)) << std::endl;
  std::cout << " " << std::endl;

  // Body to Inertial
  Mat6d b2i_spice = spice::GetFrameConversionMat(t_tdb, fixed_frame_str, "J2000");
  std::cout << "PLANET FIXED to GCRF (SPICE)" << std::endl;
  // Print with clean formatting
  std::cout << b2i_spice.format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  Mat6 b2i = RotPosVelBodyFixedToInertial(id, t_tdb);
  std::cout << "PLANET FIXED to GCRF (lupnt)" << std::endl;
  std::cout << b2i.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  std::cout << " " << std::endl;
  std::cout << " --------------------------------------------" << std::endl;
  std::cout << " " << std::endl;

  // Inertial to Body
  Mat6d i2b_spice = spice::GetFrameConversionMat(t_tdb, "J2000", fixed_frame_str);
  std::cout << "GCRF to PLANET FIXED (SPICE)" << std::endl;
  // Print with clean formatting
  std::cout << i2b_spice.format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  Mat6 i2b = RotPosVelInertialToBodyFixed(id, t_tdb);
  std::cout << "GCRF to PLANET FIXED (lupnt)" << std::endl;
  std::cout << i2b.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  return 0;
}
