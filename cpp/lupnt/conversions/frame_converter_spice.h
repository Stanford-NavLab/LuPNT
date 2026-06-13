#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/states/state.h"

namespace lupnt {

  class Cart6;

  namespace spice {
    // Vec = func(real, Vec)

    /// @brief Converts a Cartesian position+velocity state at epoch `t_tdb` from
    /// `frame_in` to `frame_out` using SPICE-derived rotations/ephemerides directly
    /// (`GetFrameConversionMat`, `GetBodyPosVel`), rather than LuPNT's native
    /// analytic frame-conversion chain (`lupnt::ConvertFrame` /
    /// frame_conversions.h).
    ///
    /// This is the SPICE-based reference/ground-truth path used to validate and
    /// calibrate LuPNT's analytic conversions -- e.g. `ComputeEopFromSpice` and
    /// `InitFrameConversionFromSpice` (frame_conversions.h) fit LuPNT's native
    /// GcrfToItrf/MoonCiToPa models to reproduce this function's results over a
    /// simulation window, and tests compare the two paths directly.
    ///
    /// @param t_tdb      Epoch [s, TDB since J2000]
    /// @param rv_in      Position+velocity state in `frame_in` [m, m/s]
    /// @param frame_in   Input reference frame
    /// @param frame_out  Output reference frame
    /// @return            Position+velocity state in `frame_out` [m, m/s]
    Vec6 ConvertFrameSpice(Real t_tdb, const Vec6& rv_in, Frame frame_in, Frame frame_out);
    /// @brief Vec3 (position-only) overload of ConvertFrameSpice(): velocity is set
    /// to zero internally before the conversion.
    Vec3 ConvertFrameSpice(Real t_tdb, const Vec3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(real, Mat)

    /// @brief Row-wise Cartesian-state overload of ConvertFrameSpice(): converts
    /// each row of `rv_in` at the single epoch `t_tdb` from `frame_in` to `frame_out`.
    MatX6 ConvertFrameSpice(Real t_tdb, const MatX6& rv_in, Frame frame_in, Frame frame_out);
    /// @brief Row-wise position overload of ConvertFrameSpice(): converts each row
    /// of `r_in` at the single epoch `t_tdb` from `frame_in` to `frame_out`.
    MatX3 ConvertFrameSpice(Real t_tdb, const MatX3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(Vec, Vec)

    /// @brief Multi-epoch overload of ConvertFrameSpice(): converts the single
    /// Cartesian state `rv_in` from `frame_in` to `frame_out` independently at each
    /// epoch in `t_tdb`, returning one output row per epoch.
    MatX6 ConvertFrameSpice(VecX t_tdb, const Vec6& rv_in, Frame frame_in, Frame frame_out);
    /// @brief Multi-epoch overload of ConvertFrameSpice(): converts the single
    /// position vector `r_in` from `frame_in` to `frame_out` independently at each
    /// epoch in `t_tdb`, returning one output row per epoch.
    MatX3 ConvertFrameSpice(VecX t_tdb, const Vec3& r_in, Frame frame_in, Frame frame_out);

    // Mat = func(Vec, Mat)

    /// @brief Time-tagged-rows overload of ConvertFrameSpice(): converts row `i` of
    /// `rv_in` from `frame_in` to `frame_out` at epoch `t_tdb(i)`, for each `i`
    /// (`t_tdb.size()` must equal `rv_in.rows()`).
    MatX6 ConvertFrameSpice(VecX t_tdb, const MatX6& rv_in, Frame frame_in, Frame frame_out);
    /// @brief Time-tagged-rows overload of ConvertFrameSpice(): converts row `i` of
    /// `r_in` from `frame_in` to `frame_out` at epoch `t_tdb(i)`, for each `i`
    /// (`t_tdb.size()` must equal `r_in.rows()`).
    MatX3 ConvertFrameSpice(VecX t_tdb, const MatX3& r_in, Frame frame_in, Frame frame_out);

    /// @brief Converts a typed `Cart6` Cartesian state to `frame_out` at epoch
    /// `t_tdb` using ConvertFrameSpice(), tagging the returned `Cart6` with its new
    /// frame.
    Cart6 ConvertFrameSpice(Real t_tdb, const Cart6& state_in, Frame frame_out);

  }  // namespace spice

}  // namespace lupnt
