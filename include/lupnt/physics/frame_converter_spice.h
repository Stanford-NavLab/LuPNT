/**
 * @file FrameConverter.cpp
 * @author Stanford NAV LAB
 * @brief Coordinate conversion functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <memory>

#include "cheby.h"
#include "lupnt/core/constants.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  class CartesianOrbitState;

  namespace spice {
    // Vec = func(real, Vec)
    Vec6 ConvertFrameSpice(Real t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out);
    Vec3 ConvertFrameSpice(Real t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(real, Mat)
    Mat<-1, 6> ConvertFrameSpice(Real t_tai, const Mat<-1, 6>& rv_in, Frame frame_in,
                                 Frame frame_out);
    Mat<-1, 3> ConvertFrameSpice(Real t_tai, const Mat<-1, 3>& r_in, Frame frame_in,
                                 Frame frame_out);

    // Mat = func(Vec, Vec)
    Mat<-1, 6> ConvertFrameSpice(VecX t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out);
    Mat<-1, 3> ConvertFrameSpice(VecX t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(Vec, Mat)
    Mat<-1, 6> ConvertFrameSpice(VecX t_tai, const Mat<-1, 6>& rv_in, Frame frame_in,
                                 Frame frame_out);
    Mat<-1, 3> ConvertFrameSpice(VecX t_tai, const Mat<-1, 3>& r_in, Frame frame_in,
                                 Frame frame_out);

    CartesianOrbitState ConvertFrameSpice(Real t_tai, const CartesianOrbitState& state_in,
                                          Frame frame_out);

  }  // namespace spice

}  // namespace lupnt
