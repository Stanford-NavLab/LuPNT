/**
 * @file cheby.c
 * @author Mark Adler
 * @brief    This illustrates how the coefficients may be extracted once using
 the SPICE library, and then used after that without the SPICE library.
   Extracted from
 * @version 0.1
 * @date 2015-08-15
 * @copyright Mark Adler (c) 2023
 */

#pragma once

#include <cspice/SpiceUsr.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lupnt/core/constants.h"

namespace lupnt {
  /// @brief Evaluate a Chebyshev polynomial (and its derivative) at `x`,
  /// double-precision implementation.
  ///
  /// Evaluates the degree-(num-1) Chebyshev series with coefficients
  /// `coeff[0..num-1]` (T_n basis, constant term first), after rescaling `x`
  /// to `[-1, 1]` via the segment's midpoint/radius `scale = {scale[0],
  /// scale[1]}` (valid for `x` in `[scale[0]-scale[1], scale[0]+scale[1]]`).
  /// Used internally by `cheby_posvel` to evaluate one (x, y, or z) component
  /// of an SPK Chebyshev record.
  ///
  /// @param x     Evaluation point (e.g. ephemeris time [s past J2000])
  /// @param scale `{midpoint, radius}` of the valid interval for `x`
  /// @param coeff Chebyshev coefficients, constant term first, length `num`
  /// @param num   Number of coefficients (polynomial degree + 1)
  /// @param f     Output: polynomial value at `x`
  /// @param df    Output: derivative of the polynomial at `x` w.r.t. `x`
  void cheby_eval(double x, double* scale, double* coeff, long num, double* f, double* df);

  /// @brief Autodiff (`Real`) counterpart of `cheby_eval`.
  ///
  /// Used by `ChebyshevFitModel::Eval` (lupnt/numerics/cheby_fit.h) to
  /// evaluate fitted Chebyshev models -- e.g. SPICE-fitted Earth-orientation
  /// parameters / lunar libration angles (`InitFrameConversionFromSpice`) and
  /// `GnssConstellation`'s fitted ECI ephemeris -- while preserving exact
  /// analytic time derivatives through `x`'s autodiff seed.
  ///
  /// @param x     Evaluation point (e.g. epoch [s past J2000])
  /// @param scale `{midpoint, radius}` of the valid interval for `x`
  /// @param coeff Chebyshev coefficients, constant term first, length `num`
  /// @param num   Number of coefficients (polynomial degree + 1)
  /// @return      `(value, derivative)` of the polynomial at `x`
  Vec2 cheby_eval_ad(Real x, double* scale, double* coeff, long num);

  /// @brief Evaluate a planetary/lunar position and velocity from a raw SPK
  /// Chebyshev-position-only data segment at time `t`, double-precision
  /// implementation.
  ///
  /// Locates the data record covering `t` within `seg` (laid out per the DAF/SPK
  /// Type-2 format described in spice_cheby.cc), then calls `cheby_eval` for
  /// each of the x/y/z components (and their derivatives, for velocity).
  ///
  /// @param t   Ephemeris time [s past J2000]
  /// @param seg Raw SPK Chebyshev segment data (from `spk_extract`/`cheby_segment`)
  /// @param len Length of `seg` in doubles
  /// @param pos Output: position `[x, y, z]` [km]
  /// @param vel Output: velocity `[vx, vy, vz]` [km/s]
  /// @return    0 on success, 1 if `t` is not covered by `seg`
  int cheby_posvel(double t, double* seg, long len, double pos[3], double vel[3]);

  /// @brief Autodiff (`Real`) counterpart of `cheby_posvel`.
  ///
  /// Called by `spice::GetBodyPosVelBase` (spice.cc) -- the core of
  /// `spice::GetBodyPosVel` -- to evaluate cached `de440.bsp` Chebyshev
  /// segments (extracted once by `spk_extract` during `LoadSpiceKernel`) at
  /// an arbitrary epoch `t`, preserving `t`'s autodiff derivative in the
  /// returned velocity.
  ///
  /// @param t   Ephemeris time [s past J2000]
  /// @param seg Raw SPK Chebyshev segment data
  /// @param len Length of `seg` in doubles
  /// @return    `[r; v]` position [km] / velocity [km/s], or `Vec6::Zero()`
  ///            if `t` is outside the segment's covered time range
  Vec6 cheby_posvel_ad(Real t, double* seg, long len);

  /// @brief Verify that a raw SPK Chebyshev segment has the expected
  /// uniform-record structure before it is used by `cheby_posvel`/`cheby_posvel_ad`.
  ///
  /// Called by `cheby_segment` immediately after reading each segment from a
  /// DAF/SPK file, to guard against segfaults on malformed data.
  ///
  /// @param seg Raw SPK Chebyshev segment data
  /// @param len Length of `seg` in doubles
  /// @return    0 if `seg` is a valid uniform set of coefficient records, 1
  ///            otherwise
  int cheby_verify(double* seg, long len);

  /// @brief Descriptor for one extracted SPK Chebyshev-position-only segment
  /// (target/center/frame codes plus the raw coefficient data), as produced
  /// by `spk_extract` / `cheby_segment` and consumed by
  /// `spice::GetBodyPosVelBase` via `cheby_posvel_ad`.
  typedef struct {
    long target;  // target body code
    long center;  // center body code
    long frame;   // frame of reference code
    long len;     // length of segment in doubles
    double* seg;  // allocated segment
  } segment_t;

  /// @brief Read one Chebyshev-position-only SPK segment from an open DAF
  /// file into `s`, verifying its structure via `cheby_verify`.
  ///
  /// Called by `spk_extract` for each Type-2 (Chebyshev position-only)
  /// segment found while scanning an SPK file.
  ///
  /// @param daf DAF file handle (from `dafopr_c`)
  /// @param dc  Segment descriptor's double-precision components (start/end
  ///            epoch [s past J2000])
  /// @param ic  Segment descriptor's integer components (target/center/frame
  ///            codes, representation type, data address range)
  /// @param s   Output: populated segment descriptor (`s->seg` is heap-allocated)
  /// @return    0 on success; 1 on error (with `s->seg` left as `nullptr`)
  int cheby_segment(SpiceInt daf, SpiceDouble* dc, SpiceInt* ic, segment_t* s);

  /// @brief Extract all Chebyshev-position-only (Type 2) segments from an SPK
  /// file into a heap-allocated array of `segment_t`, without requiring the
  /// SPICE library afterwards.
  ///
  /// Called once by `spice::LoadSpiceKernel` on `de440.bsp` (planetary/lunar
  /// ephemeris), producing the segment table that `spice::GetBodyPosVelBase`
  /// searches (by target/center body codes) and evaluates via
  /// `cheby_posvel_ad` for every `spice::GetBodyPosVel` call.
  ///
  /// @param path SPK (`.bsp`) file path
  /// @param segs Output: number of segments in the returned array
  /// @return     Heap-allocated array of `*segs` segment descriptors, or
  ///             `nullptr` on error (e.g. not a valid SPK / no Chebyshev
  ///             position-only segments)
  segment_t* spk_extract(char const* path, long* segs);

  /// @brief Free the segment data allocated by `spk_extract`.
  ///
  /// @param s Array of segment descriptors returned by `spk_extract`
  /// @param n Number of segments in `s`
  void spk_free(segment_t* s, long n);

}  // namespace lupnt
