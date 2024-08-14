#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/time_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <vector>

#include "utils.cc"

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

// Utility function to create a test matrix for Mat<-1, 6>
Mat<-1, 6> CreateTestMatrix6() {
  Mat<-1, 6> rv(2, 6);
  rv.row(0) = CreateTestVector6();
  rv.row(1) = CreateTestVector6() * 2.5;
  return rv;
}

// Utility function to create a test matrix for Mat<-1, 3>
Mat<-1, 3> CreateTestMatrix3() {
  Mat<-1, 3> r(2, 3);
  r.row(0) = CreateTestVector3();
  r.row(1) = CreateTestVector3() * 2.5;
  return r;
}

// Test cases for frame conversion
TEST_CASE("Frame_Converter") {
  Vec6 rv_init = CreateTestVector6();
  Vec3 r_init = CreateTestVector3();

  Mat<-1, 6> rv_init_mat = CreateTestMatrix6();
  Mat<-1, 3> r_init_mat = CreateTestMatrix3();

  Real t_tai = GregorianToTime(2020, 1, 1);
  VecX t_tai_vec(2);
  t_tai_vec << t_tai, t_tai + 60;

  Mat<-1, 6> rv_init_rep(2, 6);
  rv_init_rep << rv_init.transpose(), rv_init.transpose();
  Mat<-1, 3> r_init_rep(2, 3);
  r_init_rep << r_init.transpose(), r_init.transpose();

  // Vec6 = func(real, Vec6)
  for (auto frame_in : frame_list) {
    Vec6 rv_converted = ConvertFrame(t_tai, rv_init, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Vec6 rv_back = ConvertFrame(t_tai, rv_converted, frame_in, frame_out);
      Vec6 rv_final = ConvertFrame(t_tai, rv_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_VEC(rv_init, rv_final, epsilon);
    }
  }

  // Vec3 = func(real, Vec3)
  for (auto frame_in : frame_list) {
    Vec3 r_converted = ConvertFrame(t_tai, r_init, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Vec3 r_back = ConvertFrame(t_tai, r_converted, frame_in, frame_out);
      Vec3 r_final = ConvertFrame(t_tai, r_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_VEC(r_init, r_final, epsilon);
    }
  }

  // Mat<-1, 6> = func(real, Mat<-1, 6>)
  for (auto frame_in : frame_list) {
    Mat<-1, 6> rv_matrix_converted = ConvertFrame(t_tai, rv_init_mat, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Mat<-1, 6> rv_matrix_back = ConvertFrame(t_tai, rv_matrix_converted, frame_in, frame_out);
      Mat<-1, 6> rv_matrix_final = ConvertFrame(t_tai, rv_matrix_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_MAT(rv_init_mat, rv_matrix_final, epsilon);
    }
  }

  // Mat<-1, 3> = func(real, Mat<-1, 3>)
  for (auto frame_in : frame_list) {
    Mat<-1, 3> r_matrix_converted = ConvertFrame(t_tai, r_init_mat, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Mat<-1, 3> r_matrix_back = ConvertFrame(t_tai, r_matrix_converted, frame_in, frame_out);
      Mat<-1, 3> r_matrix_final = ConvertFrame(t_tai, r_matrix_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_MAT(r_init_mat, r_matrix_final, epsilon);
    }
  }

  // Mat<-1, 6> = func(VecX, Vec6)
  for (auto frame_in : frame_list) {
    Mat<-1, 6> rv_matrix_converted = ConvertFrame(t_tai_vec, rv_init, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Mat<-1, 6> rv_matrix_back = ConvertFrame(t_tai_vec, rv_matrix_converted, frame_in, frame_out);
      Mat<-1, 6> rv_matrix_final = ConvertFrame(t_tai_vec, rv_matrix_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_MAT(rv_init_rep, rv_matrix_final, epsilon);
    }
  }

  // Mat<-1, 3> = func(VecX, Vec3)
  for (auto frame_in : frame_list) {
    Mat<-1, 3> r_matrix_converted = ConvertFrame(t_tai_vec, r_init, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Mat<-1, 3> r_matrix_back = ConvertFrame(t_tai_vec, r_matrix_converted, frame_in, frame_out);
      Mat<-1, 3> r_matrix_final = ConvertFrame(t_tai_vec, r_matrix_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_MAT(r_init_rep, r_matrix_final, epsilon);
    }
  }

  // Mat<-1, 6> = func(VecX, Mat<-1, 6>)
  for (auto frame_in : frame_list) {
    Mat<-1, 6> rv_matrix_converted = ConvertFrame(t_tai_vec, rv_init_mat, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Mat<-1, 6> rv_matrix_back = ConvertFrame(t_tai_vec, rv_matrix_converted, frame_in, frame_out);
      Mat<-1, 6> rv_matrix_final = ConvertFrame(t_tai_vec, rv_matrix_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_MAT(rv_init_mat, rv_matrix_final, epsilon);
    }
  }

  // Mat<-1, 3> = func(VecX, Mat<-1, 3>)
  for (auto frame_in : frame_list) {
    Mat<-1, 3> r_matrix_converted = ConvertFrame(t_tai_vec, r_init_mat, Frame::ITRF, frame_in);
    for (auto frame_out : frame_list) {
      Mat<-1, 3> r_matrix_back = ConvertFrame(t_tai_vec, r_matrix_converted, frame_in, frame_out);
      Mat<-1, 3> r_matrix_final = ConvertFrame(t_tai_vec, r_matrix_back, frame_out, Frame::ITRF);
      REQUIRE_NEAR_REAL_MAT(r_init_mat, r_matrix_final, epsilon);
    }
  }
}
