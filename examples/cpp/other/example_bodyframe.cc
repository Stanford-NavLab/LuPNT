#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  NaifId id = NaifId::MARS;
  Real t_tdb = 100.0;
  Real t_tai = ConvertTime(t_tdb, TimeSys::TDB, TimeSys::TAI);
  Vec4 angles = PlanetOrientation(id, t_tdb);

  std::cout << "MARS Orientation at t_tdb = " << t_tdb << " s" << std::endl;
  std::cout << "  alpha0 = " << DEG * angles(0) << std::endl;
  std::cout << "  delta0 = " << DEG * angles(1) << std::endl;
  std::cout << "  W = " << DEG * angles(2) << std::endl;
  std::cout << "  Wdot = " << DEG * angles(3) << std::endl;
  std::cout << " " << std::endl;

  // Body to Inertial
  Mat6d b2i_spice =
      GetFrameConversionMat(t_tai, Frame::MARS_FIXED, Frame::GCRF);
  std::cout << "MARS_FIXED to GCRF (SPICE)" << std::endl;
  // Print with clean formatting
  std::cout << b2i_spice.format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  Mat6 b2i = RotPosVelBodyFixedToInertial(id, t_tdb);
  std::cout << "MARS_FIXED to GCRF (lupnt)" << std::endl;
  std::cout << b2i.format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  std::cout << " " << std::endl;
  std::cout << " --------------------------------------------" << std::endl;
  std::cout << " " << std::endl;

  // Inertial to Body
  Mat6d i2b_spice =
      GetFrameConversionMat(t_tai, Frame::GCRF, Frame::MARS_FIXED);
  std::cout << "GCRF to MARS_FIXED (SPICE)" << std::endl;
  // Print with clean formatting
  std::cout << i2b_spice.format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  Mat6 i2b = RotPosVelInertialToBodyFixed(id, t_tdb);
  std::cout << "GCRF to MARS_FIXED (lupnt)" << std::endl;
  std::cout << i2b.format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", ""));

  return 0;
}