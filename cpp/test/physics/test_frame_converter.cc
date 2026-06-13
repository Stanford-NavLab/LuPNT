#include <lupnt/conversions/frame_converter.h>
#include <lupnt/conversions/time_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

// Define the list of frames to test
std::vector<Frame> frame_list = {
    Frame::ITRF, Frame::GCRF,    Frame::EME,     Frame::ICRF,
    Frame::EMR,  Frame::MOON_CI, Frame::MOON_PA, Frame::MOON_ME,
};

Vec6 CreateTestVector6() { return Vec6(7000.0, 0.0, 1300.0, 0.0, 7.8, 0.0); }

Vec3 CreateTestVector3() { return Vec3(7000.0, 1300.0, 0.0); }

// Utility function to create a test matrix for MatX6
MatX6 CreateTestMatrix6() {
  MatX6 rv(2, 6);
  rv.row(0) = CreateTestVector6();
  rv.row(1) = CreateTestVector6() * 2.5;
  return rv;
}

// Utility function to create a test matrix for MatX3
MatX3 CreateTestMatrix3() {
  MatX3 r(2, 3);
  r.row(0) = CreateTestVector3();
  r.row(1) = CreateTestVector3() * 2.5;
  return r;
}

// Test cases for frame conversion
TEST_CASE("physics.frame_converter") {
  Vec6 rv_init = CreateTestVector6();
  Vec3 r_init = CreateTestVector3();

  MatX6 rv_init_mat = CreateTestMatrix6();
  MatX3 r_init_mat = CreateTestMatrix3();

  Real t_tdb = ConvertTime(GregorianToTime(2020, 1, 1), Time::TAI, Time::TDB);
  VecX t_tdb_vec(2);
  t_tdb_vec << t_tdb, t_tdb + 60;

  MatX6 rv_init_rep(2, 6);
  rv_init_rep << rv_init.transpose(), rv_init.transpose();
  MatX3 r_init_rep(2, 3);
  r_init_rep << r_init.transpose(), r_init.transpose();

  // Vec6 = func(real, Vec6)
  for (auto frame_in : frame_list) {
    Vec6 rv_in = ConvertFrame(t_tdb, rv_init, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      Vec6 rv_out = ConvertFrame(t_tdb, rv_in, frame_in, frame_out);
      Vec6 rv_final = ConvertFrame(t_tdb, rv_out, frame_out, Frame::GCRF);
      RequireNear(rv_init, rv_final, epsilon);
    }
  }

  // Vec3 = func(real, Vec3)
  for (auto frame_in : frame_list) {
    Vec3 r_in = ConvertFrame(t_tdb, r_init, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      Vec3 r_out = ConvertFrame(t_tdb, r_in, frame_in, frame_out);
      Vec3 r_final = ConvertFrame(t_tdb, r_out, frame_out, Frame::GCRF);
      RequireNear(r_init, r_final, epsilon);
    }
  }

  // MatX6 = func(real, MatX6)
  for (auto frame_in : frame_list) {
    MatX6 rv_matrix_in = ConvertFrame(t_tdb, rv_init_mat, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      MatX6 rv_matrix_out = ConvertFrame(t_tdb, rv_matrix_in, frame_in, frame_out);
      MatX6 rv_matrix_final = ConvertFrame(t_tdb, rv_matrix_out, frame_out, Frame::GCRF);
      RequireNear(rv_init_mat, rv_matrix_final, epsilon);
    }
  }

  // MatX3 = func(real, MatX3)
  for (auto frame_in : frame_list) {
    MatX3 r_matrix_in = ConvertFrame(t_tdb, r_init_mat, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      MatX3 r_matrix_out = ConvertFrame(t_tdb, r_matrix_in, frame_in, frame_out);
      MatX3 r_matrix_final = ConvertFrame(t_tdb, r_matrix_out, frame_out, Frame::GCRF);
      RequireNear(r_init_mat, r_matrix_final, epsilon);
    }
  }

  // MatX6 = func(VecX, Vec6)
  for (auto frame_in : frame_list) {
    MatX6 rv_matrix_in = ConvertFrame(t_tdb_vec, rv_init, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      MatX6 rv_matrix_out = ConvertFrame(t_tdb_vec, rv_matrix_in, frame_in, frame_out);
      MatX6 rv_matrix_final = ConvertFrame(t_tdb_vec, rv_matrix_out, frame_out, Frame::GCRF);
      RequireNear(rv_init_rep, rv_matrix_final, epsilon);
    }
  }

  // MatX3 = func(VecX, Vec3)
  for (auto frame_in : frame_list) {
    MatX3 r_matrix_in = ConvertFrame(t_tdb_vec, r_init, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      MatX3 r_matrix_out = ConvertFrame(t_tdb_vec, r_matrix_in, frame_in, frame_out);
      MatX3 r_matrix_final = ConvertFrame(t_tdb_vec, r_matrix_out, frame_out, Frame::GCRF);
      RequireNear(r_init_rep, r_matrix_final, epsilon);
    }
  }

  // MatX6 = func(VecX, MatX6)
  for (auto frame_in : frame_list) {
    MatX6 rv_matrix_in = ConvertFrame(t_tdb_vec, rv_init_mat, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      MatX6 rv_matrix_out = ConvertFrame(t_tdb_vec, rv_matrix_in, frame_in, frame_out);
      MatX6 rv_matrix_final = ConvertFrame(t_tdb_vec, rv_matrix_out, frame_out, Frame::GCRF);
      RequireNear(rv_init_mat, rv_matrix_final, epsilon);
    }
  }

  // MatX3 = func(VecX, MatX3)
  for (auto frame_in : frame_list) {
    MatX3 r_matrix_in = ConvertFrame(t_tdb_vec, r_init_mat, Frame::GCRF, frame_in);
    for (auto frame_out : frame_list) {
      MatX3 r_matrix_out = ConvertFrame(t_tdb_vec, r_matrix_in, frame_in, frame_out);
      MatX3 r_matrix_final = ConvertFrame(t_tdb_vec, r_matrix_out, frame_out, Frame::GCRF);
      RequireNear(r_init_mat, r_matrix_final, epsilon);
    }
  }
}
