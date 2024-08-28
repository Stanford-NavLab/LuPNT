#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  std::string body = "Mars";
  NaifId id = NaifId::MARS;
  Frame fixed_frame = Frame::MARS_FIXED;

  //   std::string body = "Venus";
  //   NaifId id = NaifId::VENUS;
  //   Frame fixed_frame = Frame::VENUS_FIXED;

  Real t_tdb = 20 * DAYS_YEAR * SECS_DAY;
  Real t_tai = ConvertTime(t_tdb, Time::TDB, Time::TAI);
  Vec4 angles = PlanetOrientation(id, t_tdb);

  std::cout << "Planet: " << body << std::endl;
  std::cout << "Epoch:  "
            << "TDB=" << t_tdb << " TAI=" << t_tai << std::endl;

  std::cout << "Orientation at t_tdb = " << t_tdb << " s" << std::endl;
  std::cout << "  alpha0 = " << DEG * angles(0) << std::endl;
  std::cout << "  delta0 = " << DEG * angles(1) << std::endl;
  std::cout << "  W = " << DEG * angles(2) << std::endl;
  std::cout << "  Wdot = " << DEG * angles(3) << std::endl;
  std::cout << " " << std::endl;

  // Body to Inertial
  Mat6d b2i_spice = spice::GetFrameConversionMat(t_tai, fixed_frame, Frame::GCRF);
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
  Mat6d i2b_spice = spice::GetFrameConversionMat(t_tai, Frame::GCRF, fixed_frame);
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
