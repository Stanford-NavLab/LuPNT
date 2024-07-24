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

#include "lupnt/physics/cheby.h"

#include "lupnt/numerics/math_utils.h"

namespace lupnt {
  /*
      DAF/SPK file format notes:

      ic[0] - target code
      ic[1] - center code
      ic[2] - frame code
      ic[3] - representation code (2 == Chebyshev position only)
      ic[4] - initial address of array
      ic[5] - final address of array

      len = ic[5] - ic[4] + 1

      dc[0] - initial epoch of data (seconds relative to J2000)
      dc[1] - final epoch of data (seconds relative to J2000)

      seg[len-4] - initial epoch of the first record (seconds relative to J2000)
      seg[len-3] - interval length of each record (seconds)
      seg[len-2] - elements in each record
      seg[len-1] - number of records

      seg[len-1] * seg[len-2] + 4 == len
      n = seg[len-2]
      num = (n - 2) / 3

      rec[0] - midpoint of interval covered by record (seconds relative to J2000)
      rec[1] - radius of interval (seconds)
      rec[2..num+1] - X coefficients (constant term first)
      rec[num+2..2*num+1] - Y coefficients
      rec[2*num+2..n-1] - Z coefficients

      For t, evaluate the Chebyshev polynomials T_n at (t - rec[0]) / rec[1],
      multiply by the coefficients, and sum.  The derivatives of the polynomials
      can be used to compute the velocity.  See cheby_eval() here.  The results
      are in km and km/s.  Note that all times are ephemeris times, and so do not
      take into account leap seconds.

      A SpiceDouble is simply a C double.  A SpiceInt is an integer type whose
      size is half that of double, so that two SpiceInt's fit in a SpiceDouble.
   */

  /**
   * @brief Evaluate the given Chebyshev polynomial at x, returning both the
   evaluated polynomial in *f, and the evaluated derivative of the polymonial in
   *df. The number of coefficients is num (the degree of the polynomial is num -
   1), and the coefficients are coeff[0..num-1].  The first coefficient coeff[0]
   is the constant term.  The scaling of x is provided by the midpoint scale[0]
   and the radius scale[1].  x must fall in the range scale[0] - scale[1] to
     scale[0] + scale[1].  Outside of that range, the polynomial is not valid.
  */
  void cheby_eval(double x, double* scale, double* coeff, long num, double* f, double* df) {
    double x2, w0 = 0., w1 = 0., dw0 = 0., dw1 = 0., tmp;

    x = (x - scale[0]) / scale[1];
    x2 = x * 2.;
    while (--num) {
      tmp = dw1;
      dw1 = dw0;
      dw0 = w0 * 2. + dw0 * x2 - tmp;
      tmp = w1;
      w1 = w0;
      w0 = coeff[num] + (x2 * w0 - tmp);
    }
    *f = coeff[0] + (x * w0 - w1);
    *df = (w0 + x * dw0 - dw1) / scale[1];
  }

  Vec2 cheby_eval_ad(Real x, double* scale, double* coeff, long num) {
    Real x2, w0 = 0., w1 = 0., dw0 = 0., dw1 = 0., tmp;
    x = (x - scale[0]) / scale[1];
    x2 = x * 2.;
    while (--num) {
      tmp = dw1;
      dw1 = dw0;
      dw0 = w0 * 2. + dw0 * x2 - tmp;
      tmp = w1;
      w1 = w0;
      w0 = coeff[num] + (x2 * w0 - tmp);
    }

    Vec2 ret_state{coeff[0] + (x * w0 - w1), (w0 + x * dw0 - dw1) / scale[1]};
    return ret_state;
  }

  /**
   * @brief Find the appropriate SPK record for time t and compute the position
   * and velocity for that time.  Returns 0 on success, 1 if the time is not
   * covered by the segment. */
  int cheby_posvel(double t, double* seg, long len, double pos[3], double vel[3]) {
    long k, num;

    k = (long)floor((t - seg[len - 4]) /   // seg[len-4] is initial epoch
                    seg[len - 3]);         // seg[len-3] is record span
    if (k < 0 || k >= (long)seg[len - 1])  // seg[len-1] is number of records
      return 1;
    num = (long)seg[len - 2];  // seg[len-2] is size of record
    seg += k * num;            // point seg to the record for t
    num = (num - 2) / 3;       // number of coefficients
    cheby_eval(t, seg, seg + 2, num, pos, vel);
    cheby_eval(t, seg, seg + 2 + num, num, pos + 1, vel + 1);
    cheby_eval(t, seg, seg + 2 + 2 * num, num, pos + 2, vel + 2);
    return 0;
  }

  Vec6 cheby_posvel_ad(Real t, double* seg, long len) {
    long k, num;

    k = (long)floor((t.val() - seg[len - 4]) /  // seg[len-4] is initial epoch
                    seg[len - 3]);              // seg[len-3] is record span
    if (k < 0 || k >= (long)seg[len - 1])       // seg[len-1] is number of records
      return Vec6::Zero();

    num = (long)seg[len - 2];  // seg[len-2] is size of record
    seg += k * num;            // point seg to the record for t
    num = (num - 2) / 3;       // number of coefficients

    Vec2 xdx = cheby_eval_ad(t, seg, seg + 2, num);
    Vec2 ydy = cheby_eval_ad(t, seg, seg + 2 + num, num);
    Vec2 zdz = cheby_eval_ad(t, seg, seg + 2 + 2 * num, num);

    Vec6 posvel{xdx[0], ydy[0], zdz[0], xdx[1], ydy[1], zdz[1]};
    return posvel;
  }

  /**
   * @brief Verify that the provided segment meets the constraints of a uniform
   set of coefficient records.  Return 0 on success or 1 if the segment is
   invalid. This should be done before using the segment in order to avoid
   segfaults on invalid data. */
  int cheby_verify(double* seg, long len) {
    double recs = seg[len - 1],  // number of records
        elts = seg[len - 2],     // elements (doubles) in each record
        span = seg[len - 3],     // time span of each record in seconds
        init = seg[len - 4];     // initial epoch in seconds relative to J2000
    long n, k;
    double *p, *q;

    if (recs != (long)recs ||                      // recs is an integer
        elts != (long)elts ||                      // elts is an integer
        (long)recs * (long)elts + 4 != len ||      // total length is correct
        3 * (((long)elts - 2) / 3) + 2 != elts ||  // integer number of coeffs
        seg[0] - seg[1] != init ||                 // 1st start is init
        span != 2 * seg[1])                        // 1st radius matches span
      return 1;
    n = (long)recs;
    k = (long)elts;
    p = seg;
    while (--n) {
      q = p + k;                       // scan all q following p
      if (q[1] != p[1] ||              // all radii the same
          q[0] - q[1] != p[0] + p[1])  // next start is last end
        return 1;
      p = q;
    }
    return 0;
  }

  /**
   * @brief Print an error message
   */
  void cheby_err(char const* msg, ...) {
    fputs("cheby error: ", stderr);
    va_list ap;
    va_start(ap, msg);
    vfprintf(stderr, msg, ap);
    va_end(ap);
    putc('\n', stderr);
  }

  /**
   * @brief Load one segment of an SPK file, which covers one target over a range
   of epochs.  Save the target code, reference location code for the target
     position, and the reference frame code.  Load the segment and verify its
     structure.  On success return 0.  If there is an error, return 1 and set
     s->seg to NULL. */
  int cheby_segment(SpiceInt daf, SpiceDouble* dc, SpiceInt* ic, segment_t* s) {
    SpiceDouble* last;

    // save segment codes
    s->target = ic[0];
    s->center = ic[1];
    s->frame = ic[2];

    // allocate memory for the segment and read it in
    s->len = ic[5] - ic[4] + 1;  // number of doubles in segment
    // s->seg = malloc(s->len * sizeof(SpiceDouble));
    s->seg = (double*)malloc(s->len * sizeof(double));

    if (s->seg == NULL) {
      cheby_err("out of memory");
      return 1;
    }
    dafgda_c(daf, ic[4], ic[5], s->seg);  // load segment
    if (failed_c()) {
      reset_c();
      free(s->seg);
      s->seg = NULL;
      cheby_err("could not read SPK segment from file");
      return 1;
    }

    // verify the integrity of the segment
    last = s->seg + s->len - 4 - (long)(s->seg[s->len - 2]);
    if (cheby_verify(s->seg, s->len) ||  // segment structure ok
        dc[0] != s->seg[s->len - 4] ||   // start epoch matches
        dc[1] != last[0] + last[1]) {    // end epoch matches
      free(s->seg);
      s->seg = NULL;
      cheby_err("SPK segment format is invalid");
      return 1;
    }

    // return loaded segment
    return 0;
  }

  /**
   * @brief Scan through the SPK file path and extract all of the Chebyshev
     position-only segments, saving them in an allocated array of segment_t,
     which is returned.  If there is an error, NULL is returned.  *segs is set to
     the number of segments in the array.  Once this is done, this array can be
     used by cheby_verify() and cheby_posvel() above, with no dependency on or
     reference to the SPICE library.
  */
  segment_t* spk_extract(char const* path, long* segs) {
    SpiceInt daf;
    SpiceBoolean found;
    union {
      SpiceDouble d[128];
      SpiceChar c[1024];
    } sum;
    const SpiceInt nd = 2, ni = 6;
    SpiceDouble dc[nd];
    SpiceInt ic[ni];
    segment_t *spk, *mem;

    // turn off error reporting and aborts for SPICE functions
    errprt_c("set", 0, (char*)"none");
    erract_c("set", 0, (char*)"return");

    // open the file and verifiy that it is a DAF SPK file
    dafopr_c(path, &daf);  // open SPK file for reading
    if (failed_c()) {
      reset_c();
      cheby_err("could not open %s as a DAF", path);
      return NULL;
    }
    dafgsr_c(daf, 1, 1, 128, sum.d, &found);  // read first record
    if (failed_c() || !found || memcmp(sum.c, "DAF/SPK ", 8)) {
      reset_c();
      dafcls_c(daf);
      cheby_err("%s is not an SPK file", path);
      return NULL;
    }

    // count the number of Chebyshev position-only segments in the DAF file
    *segs = 0;
    dafbfs_c(daf);                     // begin forward search
    while (daffna_c(&found), found) {  // find the next array
      dafgs_c(sum.d);                  // get array summary
      dafus_c(sum.d, nd, ni, dc, ic);  // unpack the array summary
      if (failed_c()) break;
      if (ic[3] == 2)  // Chebyshev position only
        (*segs)++;     // count segment
    }
    if (failed_c() || *segs == 0) {
      reset_c();
      dafcls_c(daf);
      cheby_err("file error or Chebyshev position-only segments in %s", path);
      return NULL;
    }

    // allocate table of segment descriptors
    spk = (segment_t*)malloc(*segs * sizeof(segment_t));
    if (spk == NULL) {
      dafcls_c(daf);
      cheby_err("out of memory");
      return NULL;
    }

    // read and save the Chebyshev position-only segments
    *segs = 0;
    dafbfs_c(daf);                     // begin forward search
    while (daffna_c(&found), found) {  // find the next array
      dafgs_c(sum.d);                  // get array summary
      dafus_c(sum.d, nd, ni, dc, ic);  // unpack the array summary
      if (failed_c()) break;
      if (ic[3] == 2 && !cheby_segment(daf, dc, ic, spk + *segs)) (*segs)++;
    }
    if (failed_c() || *segs == 0) {
      reset_c();
      dafcls_c(daf);
      free(segs);
      cheby_err("no valid Chebyshev position-only segments in %s", path);
      return NULL;
    }

    // close the DAF file and return segment table
    dafcls_c(daf);
    errprt_c("set", 0, (char*)"short");
    erract_c("set", 0, (char*)"abort");
    mem = (segment_t*)realloc(spk, *segs * sizeof(segment_t));
    if (mem != NULL) spk = mem;
    return spk;
  }

  /**
   * @brief
   * Free the resources of an SPK structure created by spk_extract(). */
  void spk_free(segment_t* s, long n) {
    long i;

    for (i = 0; i < n; i++) free(s[i].seg);
    free(s);
  }

}  // namespace lupnt
