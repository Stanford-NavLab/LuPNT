/**
 * @file ChebyshevExample.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-02-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <lupnt/physics/SpiceInterface.h>

/**
 * @brief
 * Load the SPK files on the command line and verify the position and velocity
   at J2000 + 0 seconds for each Chebyshev position-only segment against the
   same result from the SPICE library. */

using namespace LPT;

int main() {
  char const *filepath = "../data/ephemeris/de440.bsp";

  segment_t *s;
  long n, i;
  double pos[3], vel[3];
  SpiceInt eph, frame, center;
  SpiceDouble desc[5], pv[6];
  SpiceBoolean found;
  SpiceChar id[41];

  s = spk_extract(filepath, &n);
  if (s == NULL) {
    cheby_err("could not load %s as an SPK file", filepath);
  }
  furnsh_c(filepath);
  for (i = 0; i < n; i++) {
    // show segment info and position and velocity at J2000 + 0
    printf("target = %ld, center = %ld, frame = %ld\n", s[i].target,
           s[i].center, s[i].frame);
    if (s[i].seg == NULL || cheby_verify(s->seg, s->len)) {
      cheby_err("bad segment");
      putchar('\n');
    }
    if (cheby_posvel(0, s[i].seg, s[i].len, pos, vel)) {
      cheby_err("J2000 + 0 out of range (!)");
      putchar('\n');
    }
    printf("pos(0) = (%g, %g, %g)\n", pos[0], pos[1], pos[2]);
    printf("vel(0) = (%g, %g, %g)\n", vel[0], vel[1], vel[2]);

    // check position and velocity against SPICE library access
    spksfs_c(s[i].target, 0, sizeof(id), &eph, desc, id, &found);
    if (!found) {
      cheby_err("target %d not found!", s[i].target);
      putchar('\n');
    }
    spkpvn_c(eph, desc, 0, &frame, pv, &center);
    if (s[i].frame != frame || s[i].center != center)
      cheby_err("codes mismatch");
    if (pos[0] != pv[0] || pos[1] != pv[1] || pos[2] != pv[2])
      cheby_err("position mismatch");
    if (vel[0] != pv[3] || vel[1] != pv[4] || vel[2] != pv[5])
      cheby_err("velocity mismatch");
    putchar('\n');
  }
  unload_c(filepath);
  spk_free(s, n);
  return 0;
}