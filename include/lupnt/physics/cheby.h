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
  /**
  * @brief Evaluate the given Chebyshev polynomial at x, returning both the
  evaluated polynomial in *f, and the evaluated derivative of the polymonial in
  *df. The number of coefficients is num (the degree of the polynomial is num -
  1), and the coefficients are coeff[0..num-1].  The first coefficient coeff[0] is
  the constant term.  The scaling of x is provided by the midpoint scale[0] and
  the radius scale[1].  x must fall in the range scale[0] - scale[1] to
  scale[0] + scale[1].  Outside of that range, the polynomial is not valid.
  */
  void cheby_eval(double x, double* scale, double* coeff, long num, double* f, double* df);
  Vec2 cheby_eval_ad(Real x, double* scale, double* coeff, long num);

  /**
   * @brief Find the appropriate SPK record for time t and compute the position
   * and velocity for that time.  Returns 0 on success, 1 if the time is not
   * covered by the segment. */
  int cheby_posvel(double t, double* seg, long len, double pos[3], double vel[3]);
  Vec6 cheby_posvel_ad(Real t, double* seg, long len);

  /**
   * @brief Verify that the provided segment meets the constraints of a uniform
   set of coefficient records.  Return 0 on success or 1 if the segment is
   invalid. This should be done before using the segment in order to avoid
   segfaults on invalid data. */
  int cheby_verify(double* seg, long len);

  /**
   * @brief SPK segment descriptor.
   * */
  typedef struct {
    long target;  // target body code
    long center;  // center body code
    long frame;   // frame of reference code
    long len;     // length of segment in doubles
    double* seg;  // allocated segment
  } segment_t;

  /**
   * @brief Load one segment of an SPK file, which covers one target over a range
   of epochs.  Save the target code, reference location code for the target
     position, and the reference frame code.  Load the segment and verify its
     structure.  On success return 0.  If there is an error, return 1 and set
     s->seg to NULL. */
  int cheby_segment(SpiceInt daf, SpiceDouble* dc, SpiceInt* ic, segment_t* s);

  /**
   * @brief Scan through the SPK file path and extract all of the Chebyshev
     position-only segments, saving them in an allocated array of segment_t,
     which is returned.  If there is an error, NULL is returned.  *segs is set to
     the number of segments in the array.  Once this is done, this array can be
     used by cheby_verify() and cheby_posvel() above, with no dependency on or
     reference to the SPICE library.
  */
  segment_t* spk_extract(char const* path, long* segs);

  /**
   * @brief
   * Free the resources of an SPK structure created by spk_extract(). */
  void spk_free(segment_t* s, long n);

}  // namespace lupnt
