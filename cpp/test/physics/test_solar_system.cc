#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../data.cc"
#include "../utils.cc"
using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("Solar System") {
  // Time
  Real t_tai = GregorianToTime("2025-01-01T12:00:00.000");
  Real t_tt = ConvertTime(t_tai, Time::TAI, Time::TT);
  Real t_ut1 = ConvertTime(t_tai, Time::TAI, Time::UT1);
  Real t_tdb = ConvertTime(t_tai, Time::TAI, Time::TDB);
  Real mjd_tt = TimeToMjd(t_tt);
  Real mjd_ut1 = TimeToMjd(t_ut1);

  // Mean obliquity
  Real eps0 = MeanObliquity(mjd_tt);

  // Equatorial to ecliptic matrix
  Mat3 R_eq2ecl = Equatorial2EclipticMatrix(mjd_tt);

  // Precession matrix
  Mat3 R_pre = PrecessionMatrix(MJD_J2000_TT, mjd_tt);

  // Nutation angles
  auto [dpsi, deps] = NutAngles(mjd_tt);
  // Nutation matrix
  Mat3 R_nut = NutationMatrix(mjd_tt);
  // Nutation matrix low precision
  Mat3 R_nut_low = NutationMatrixLowPrecision(mjd_tt);
  // Equation of the equinoxes
  Real eqeq = EquinoxEquation(mjd_tt);
  // Sun position
  Vec3 r_sun = SunPositionLowPrecision(mjd_tt);
  // Moon position
  Vec3 r_moon = MoonPositionLowPrecision(mjd_tt);
  // Greenwich hour angle matrix
  Mat3 R_gha = GreenwichHourAngleMatrix(mjd_ut1);

  // Planet orientation
  Vec4 angles_mercury = PlanetOrientation(BodyId::MERCURY, t_tdb);
  Vec4 angles_venus = PlanetOrientation(BodyId::VENUS, t_tdb);
  Vec4 angles_mars = PlanetOrientation(BodyId::MARS, t_tdb);
  Vec4 angles_jupiter = PlanetOrientation(BodyId::JUPITER, t_tdb);
  Vec4 angles_saturn = PlanetOrientation(BodyId::SATURN, t_tdb);
  Vec4 angles_uranus = PlanetOrientation(BodyId::URANUS, t_tdb);
  Vec4 angles_neptune = PlanetOrientation(BodyId::NEPTUNE, t_tdb);

  // Rotation matrix from inertial to body-fixed frame
  Mat3 R_i2b_mercury = RotPosInertialToBodyFixed(BodyId::MERCURY, t_tdb);
  Mat3 R_i2b_venus = RotPosInertialToBodyFixed(BodyId::VENUS, t_tdb);
  Mat3 R_i2b_mars = RotPosInertialToBodyFixed(BodyId::MARS, t_tdb);
  Mat3 R_i2b_jupiter = RotPosInertialToBodyFixed(BodyId::JUPITER, t_tdb);
  Mat3 R_i2b_saturn = RotPosInertialToBodyFixed(BodyId::SATURN, t_tdb);
  Mat3 R_i2b_uranus = RotPosInertialToBodyFixed(BodyId::URANUS, t_tdb);
  Mat3 R_i2b_neptune = RotPosInertialToBodyFixed(BodyId::NEPTUNE, t_tdb);

  // Rotation matrix from body-fixed to inertial frame
  Mat3 R_b2i_mercury = RotPosBodyFixedToInertial(BodyId::MERCURY, t_tdb);
  Mat3 R_b2i_venus = RotPosBodyFixedToInertial(BodyId::VENUS, t_tdb);
  Mat3 R_b2i_mars = RotPosBodyFixedToInertial(BodyId::MARS, t_tdb);
  Mat3 R_b2i_jupiter = RotPosBodyFixedToInertial(BodyId::JUPITER, t_tdb);
  Mat3 R_b2i_saturn = RotPosBodyFixedToInertial(BodyId::SATURN, t_tdb);
  Mat3 R_b2i_uranus = RotPosBodyFixedToInertial(BodyId::URANUS, t_tdb);
  Mat3 R_b2i_neptune = RotPosBodyFixedToInertial(BodyId::NEPTUNE, t_tdb);

  // Rotation matrix from inertial to body-fixed frame (position and velocity)
  Mat6 R_i2b_vel_mercury = RotPosVelInertialToBodyFixed(BodyId::MERCURY, t_tdb);
  Mat6 R_i2b_vel_venus = RotPosVelInertialToBodyFixed(BodyId::VENUS, t_tdb);
  Mat6 R_i2b_vel_mars = RotPosVelInertialToBodyFixed(BodyId::MARS, t_tdb);
  Mat6 R_i2b_vel_jupiter = RotPosVelInertialToBodyFixed(BodyId::JUPITER, t_tdb);
  Mat6 R_i2b_vel_saturn = RotPosVelInertialToBodyFixed(BodyId::SATURN, t_tdb);
  Mat6 R_i2b_vel_uranus = RotPosVelInertialToBodyFixed(BodyId::URANUS, t_tdb);
  Mat6 R_i2b_vel_neptune = RotPosVelInertialToBodyFixed(BodyId::NEPTUNE, t_tdb);

  // Rotation matrix from body-fixed to inertial frame (position and velocity)
  Mat6 R_b2i_vel_mercury = RotPosVelBodyFixedToInertial(BodyId::MERCURY, t_tdb);
  Mat6 R_b2i_vel_venus = RotPosVelBodyFixedToInertial(BodyId::VENUS, t_tdb);
  Mat6 R_b2i_vel_mars = RotPosVelBodyFixedToInertial(BodyId::MARS, t_tdb);
  Mat6 R_b2i_vel_jupiter = RotPosVelBodyFixedToInertial(BodyId::JUPITER, t_tdb);
  Mat6 R_b2i_vel_saturn = RotPosVelBodyFixedToInertial(BodyId::SATURN, t_tdb);
  Mat6 R_b2i_vel_uranus = RotPosVelBodyFixedToInertial(BodyId::URANUS, t_tdb);
  Mat6 R_b2i_vel_neptune = RotPosVelBodyFixedToInertial(BodyId::NEPTUNE, t_tdb);
}
