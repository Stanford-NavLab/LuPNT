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
    MatX6 ConvertFrameSpice(Real t_tai, const MatX6& rv_in, Frame frame_in, Frame frame_out);
    MatX3 ConvertFrameSpice(Real t_tai, const MatX3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(Vec, Vec)
    MatX6 ConvertFrameSpice(VecX t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out);
    MatX3 ConvertFrameSpice(VecX t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(Vec, Mat)
    MatX6 ConvertFrameSpice(VecX t_tai, const MatX6& rv_in, Frame frame_in, Frame frame_out);
    MatX3 ConvertFrameSpice(VecX t_tai, const MatX3& r_in, Frame frame_in, Frame frame_out);

    CartesianOrbitState ConvertFrameSpice(Real t_tai, const CartesianOrbitState& state_in,
                                          Frame frame_out);

  }  // namespace spice

}  // namespace lupnt
