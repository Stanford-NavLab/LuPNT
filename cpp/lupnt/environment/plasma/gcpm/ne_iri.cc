/**
 * @file ne_iri.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-06
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "lupnt/environment/plasma/gcpm/ne_iri.h"

#include <array>
#include <cmath>
#include <iostream>

#include "lupnt/environment/plasma/core/math_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {

  double ne_iri_ps_trough(double r, double al, double alatr, double amlt, double akp,
                          const std::array<int, 2>& itime) {
    double val = 0.0;

    // We don't go inside the Earth
    if (r <= 1.0) return val;

    // Altitude of point of interest
    double aheight = (r - 1.0) * RE;
    double ahemisphere = sign(1.0, alatr);

    // density at the equator for given L-shell
    double eq_iri_ps_trough = ne_iri_ps_trough_eq(al, amlt, akp, itime);

#if USE_STATIC_GCPM_CPP
    static double amlt_o = -99.0, akp_o = -99.0, al_o = -99.0;
    static int itime1_o = -99, itime2_o = -99;
    static double ahemisphere_o = -99.0;
    static double transh, rf2, alpha, dno, co, switchh, switchw;
    static int istat = 0;
#else
    double amlt_o = -99.0, akp_o = -99.0, al_o = -99.0;
    int itime1_o = -99, itime2_o = -99;
    double ahemisphere_o = -99.0;
    double transh, rf2, alpha, dno, co, switchh, switchw;
    int istat = 0;
#endif

    IRIParams params;

    if (amlt != amlt_o || akp != akp_o || itime[0] != itime1_o || itime[1] != itime2_o
        || std::abs(al - al_o) > 1.0e-5 || istat < 0 || ahemisphere != ahemisphere_o) {
      // determine the height power law fit parameters that connect the
      // topside ionosphere to the equatorial density along given L-shell
      iri_ps_bridge(r, al, alatr, amlt, itime, eq_iri_ps_trough, transh, rf2, alpha, dno, co,
                    switchh, switchw, istat);
      if (istat < 0) {
        double along = (amlt - 12.0) * AMLTRAD;
        params = iri_sm(alatr, along, r, itime);
        return params.neiri;
      }
      amlt_o = amlt;
      akp_o = akp;
      al_o = al;
      itime1_o = itime[0];
      itime2_o = itime[1];
      ahemisphere_o = ahemisphere;
    }

    // compute density as given by the bridge function
    double swtchb = (aheight <= (switchh - switchw)) ? 0.0
                    : (aheight >= (switchh + switchw))
                        ? 1.0
                        : (aheight - (switchh - switchw)) / (2.0 * switchw);

    double eq_bridge
        = (dno * pow(aheight, -alpha) + co) * (1.0 - swtchb) + swtchb * eq_iri_ps_trough;

    // to call IRI or not to call IRI, that is the question...
    if (aheight < 1000.0) {
      double along = (amlt - 12.0) * AMLTRAD;
      params = iri_sm(alatr, along, r, itime);
      if (r <= rf2) {
        return params.neiri;  // outf[0][0]
      } else {
        double swtch1 = switchon(aheight, transh, 5.0);
        val = params.neiri * (1.0 - swtch1) + eq_bridge * swtch1;
        val = val * double(istat + 1) + params.neiri * double(std::abs(istat));
        return val;
      }
    } else {
      val = eq_bridge;
    }
    return val;
  }

  void iri_ps_bridge(double rr, double al, double alatr, double amlt,
                     const std::array<int, 2>& itime, double eq_iri_ps_trough, double& transh,
                     double& rf2, double& alpha, double& dno, double& co, double& switchh,
                     double& switchw, int& istat) {
    // Constants
    double tot_delh = 600.0 / RE;
    double r = 260.0 / RE + 1.0;
    double amltrad = 3.1415927 / 12.0;

    // istat must be either 0 or -1 as it is used in an equation later
    istat = 0;
    double dl = al;

    double ahemisphere = sign(1.0, alatr);

    // get height and densiy of the f2 peak
    double along = amod((amlt + 12.0), 24.0) * AMLTRAD;
    double cosrl = std::min(std::sqrt(r / al), 1.0);
    double alatrl = std::acos(cosrl) * ahemisphere;

    IRIParams params = iri_sm(alatrl, along, r, itime);

    double r1, r2;
    for (int i = 0; i < 2; i++) {
      r2 = params.hmf2_km / RE + 1.0;
      cosrl = std::min(std::sqrt(r2 / al), 1.0);
      alatrl = std::acos(cosrl) * ahemisphere;
      params = iri_sm(alatrl, along, r2, itime);
    }

    // approximate the F2 peak along the L-shell=al
    rf2 = params.hmf2_km / RE + 1.0;

    // If L-shell is at or below "r", the starting radial distance
    // for searching for the maximum negative slope above the f2 peak,
    // then the L-shell provided is exclusively an ionospheric issue
    // and we need to pass back parameters that will minimize the hassle
    // associated with the rest of the calculation for density, which
    // will necessarily exclude the bridge density anyway.
    if (rr <= rf2) {
      istat = -1;
      return;
    }

    // In an effort to reduce the cals to iri2007 the following is used
    // to approximate the point of maximum negative slope in the topside
    // ionosphere. This has been obtained from a linear fit to this location
    // (derived from the search algorithm above) as a function of returned
    // rz12 value from IRI2007. That analysis obtained this relationship:

    // Todo: Get IRI parameters
    // iri = get_iri_params(itime)
    double ro = 1.05454 + 8.62678e-5 * params.rz12;

    if (ro <= rf2) {
      ro = rf2 + 0.01;
    }

    transh = (ro - 1.0) * RE;
    double diffh = 1.0;
    double ah1 = transh - diffh;
    double ah2 = transh + diffh;
    r1 = ah1 / RE + 1.0;
    r2 = ah2 / RE + 1.0;

    cosrl = std::min(std::sqrt(ro / al), 1.0);
    alatrl = std::acos(cosrl) * ahemisphere;
    params = iri_sm(alatrl, along, ro, itime);
    double antransh = params.neiri;  // outf[0][0]

    // setup for use of densities and heights of the locations
    // on either side of the point of maximum negative slope.
    // Since only calculating ro from a fitted function, need to separately
    // determine the ionospheric densities above and below to support initial
    // calculation of the power law function.
    cosrl = std::min(std::sqrt(r1 / al), 1.0);
    alatrl = std::acos(cosrl) * ahemisphere;
    params = iri_sm(alatrl, along, r1, itime);
    double an1 = params.neiri;  // outf[0][0]

    cosrl = std::min(std::sqrt(r2 / al), 1.0);
    alatrl = std::acos(cosrl) * ahemisphere;
    params = iri_sm(alatrl, along, r2, itime);
    double an2 = params.neiri;  // outf[0][0]

    if (al <= r2) {
      istat = -1;  // no bridge required
      return;
    }

    double eqh = (al - 1.0) * RE;
    alpha = -std::log10(an1 / an2) / std::log10(ah1 / ah2);
    double ano = an1 * pow(ah1, alpha);
    double an3 = ano * pow(eqh, -alpha);

    switchh = eqh * 2.0;
    switchw = eqh / 10.0;

    if (eq_iri_ps_trough >= an3) {
      if (an2 <= eq_iri_ps_trough) {
        alpha = std::log10(antransh / eq_iri_ps_trough) / std::log10(transh / eqh);
        ano = antransh * pow(transh, alpha);
        dno = ano;
      } else {
        co = eq_iri_ps_trough - an3;
        alpha = -std::log10((an1 - co) / (an2 - co)) / std::log10(ah1 / ah2);
        ano = (an1 - co) * pow(ah1, alpha);
        dno = ano;
      }
    } else {
      switchh = transh + (eqh - transh) / 2.0;
      switchw = (eqh - transh) / 2.0;
      dno = ano;
      co = 0.0;
    }

    return;
  }

  double ne_iri_ps_trough_eq(double al, double amlt, double akp, const std::array<int, 2>& itime) {
#if USE_STATIC_GCPM_CPP
    static double amlt_o = -99.0, akp_o = -99.0;
    static int itime1_o = -99, itime2_o = -99;
    static double a8;
    static double am1, b1, x234;
    static double transh, alpha, ano, rintercept;
#else
    double amlt_o = -99.0, akp_o = -99.0;
    int itime1_o = -99, itime2_o = -99;
    double a8;
    double am1, b1, x234;
    double transh, alpha, ano, rintercept;
#endif

    double geosync_trough = 0.0;

    if (al <= 1.0) {
      return 0.0;
    }

    double aheight = (al - 1.0) * RE;

    double pp_factor = pp_profile(al, amlt, akp, a8);

    double ps_inner = ne_inner_ps(al, amlt, itime, am1, b1, x234) * 1.0e6;

    if (amlt != amlt_o || akp != akp_o || itime[0] != itime1_o || itime[1] != itime2_o) {
      iri_ps_eq_bridge(al, amlt, itime, transh, alpha, ano, am1, b1, x234, rintercept);
      amlt_o = amlt;
      akp_o = akp;
      itime1_o = itime[0];
      itime2_o = itime[1];
    }

    double ps_bridge = ano * pow(aheight, -alpha);
    double off = 0.0;
    double transwidth = 0.02;
    double swtch2 = switchon(al, rintercept + off, transwidth);
    double swtch3 = switchon(al, rintercept - off, transwidth);

    double alatr = 0.0;
    double along = (amlt - 12.0) * AMLTRAD;
    IRIParams params = iri_sm(alatr, along, al, itime);

    double swtch1 = switchon(aheight, transh, 5.0);
    double ne_eq_trough_den = ne_eq_trough(al, amlt, akp, geosync_trough);
    double zl = check_crossing(a8, am1, b1, x234, amlt, akp, geosync_trough);
    double diff = a8 - zl;
    double offset = (0.0166513 - 0.0450188 * diff) * (1.0 - switchon(diff, 0.3698744, 0.05));
    double swtch4 = switchon(al, zl + offset, 0.3);
    double swtch5 = switchon(al, zl - offset, 0.3);

    double val
        = params.neiri * (1.0 - swtch1)
          + ((ps_bridge * (1.0 - swtch2) * swtch1 + ps_inner * swtch3) * pp_factor) * (1.0 - swtch4)
          + ne_eq_trough_den * 1.0e6 * swtch5;

    return val;
  }

  double pp_profile(double al, double amlt, double akp, double& a8) {
#if USE_STATIC_GCPM_CPP
    static double akp_old = -99.0, amlt_old = -99.0;
    static double a9, centroid;
#else
    double akp_old = -99.0, amlt_old = -99.0;
    double a9, centroid;
#endif

    if ((akp != akp_old) || (amlt != amlt_old)) {
      bulge(amlt, akp, a8, a9, centroid);
    }
    akp_old = akp;
    amlt_old = amlt;

    double factor = std::min(27.75, 2.0 * (a9 - 1.0) * std::log10(al / a8));
    double pp_profile = std::pow(1.0 + std::pow(10.0, factor), -a9 / (a9 - 1.0));

    return pp_profile;
  }

  /**
   * Calculates the magnetic local time (MLT) location of the plasmapause bulge
   * centroid, as well as the parameters a8 and a9 used in modeling the
   * plasmapause transition.
   *
   * The bulge centroid calculation was corrected based on an interpretation error
   * in Moldwin et al.'s equation for the bulge centroid location in MLT.
   *
   * @param amlt Magnetic local time in hours
   * @param akp Planetary Kp index
   * @param a8 Output parameter representing the plasmapause location factor
   * @param a9 Output parameter representing the plasmapause steepness factor
   * @param centroid Output parameter representing the bulge centroid location
   */
  void bulge(double amlt, double akp, double& a8, double& a9, double& centroid) {
    const double AHOUR_RAD = 0.26179939;  // Conversion factor from hours to radians
    const double AHRRAD = 2.6179939e-1;   // Another conversion factor used in calculations

    // std::cout << "Entering bulge: " << amlt << ", " << akp << ", " << a8 << ",
    // " << a9 << ", " << centroid << std::endl;

    // Calculate the bulge centroid location based on Kp index
    centroid = 47.0 / (akp + 3.9) + 11.3;
    double x = amlt - centroid;
    if (x < -12.0) x += 24.0;
    if (x > 12.0) x -= 24.0;
    double absx = std::abs(x) * AHRRAD;

    // Compute MLT-dependent parameters using sinusoidal variations
    double along = amlt * AHOUR_RAD + 1.5707963;
    double salong = std::sin(along);
    double b1 = 0.043 * salong - 0.4589;
    double b2 = -0.361 * salong + 5.7464;

    // Compute plasmapause location factor
    a8 = (b1 * akp + b2) * (1.0 + std::exp(-1.5 * absx * absx + 0.08 * absx - 0.7));

    // Compute plasmapause steepness factor
    double b3 = -0.0243 * salong + 0.2464;
    double b4 = -0.3137 * salong - 5.2214;
    double b5 = 3.5817 * salong + 48.8114;
    a9 = b3 * akp * akp + b4 * akp + b5;
  }

  /**
   * @brief Function to compute the trough density in the equatorial plane
   * @param al (L-shell)
   * @param amlt (solar magnetic local time)
   * @param akp (Global Kp index)
   * @return double (trough density in the equatorial plane cm^-3)
   */
  double ne_eq_trough(double al, double amlt, double akp, double& geosync_trough) {
    // compute MLT where trough density peaks assuming constant Kp
    double phitp = 0.145 * akp * akp - 2.63 * akp + 21.86;

    // assuming constant density increase of 0.56+-0.08 cm-3h-1 before phitp
    // and constant density decrease of -0.83+-0.15 cm-3h-1 after that and before
    // 1 hour MLT, where the decrease is replaced by zero growth.  Growth begins
    // again at 3.5 hours MLT.  A minimum density of 0.18 cm-3 is used.
    // Transitions between growth rates are made over an hour or more, as shown
    // below. compute unsmoothed peak density at phitp
    double antp = (phitp - 3.5) * 0.56;

    // how long with it take to loose what was gained?
    //  minimum loss rate required to get back to 0.18 cm-3 by 2.0 hours MLT
    double damping_time = std::min(26.0 - phitp, antp / 0.83);

    // given the loss rate above as a minimum.
    double damping = -antp / damping_time;
    double down_time = phitp + damping_time;
    double del = 3.5 - (down_time - 24.0);
    double center = 3.5 - del / 2.0;
    if (center < 0.0) center = 24.0 + center;
    double diff = amlt - center;
    if (diff < -12.0) diff += 24.0;
    if (diff > 12.0) diff -= 24.0;
    double aminden = 0.18;
    double width = 2.0 * del;
    double denmin = aminden + diff * diff / (del * width);

    // Compute density growth or decay based on MLT
    double dengrow = 0.56 * (amlt - 3.5) + aminden;
    double sdel = 0.4, shift = 0.5;
    double switch1 = switchon(amlt, 3.5 + shift, sdel);
    double switch2 = switchon(amlt, phitp, 0.5);
    double dendamp, switch0, switch3;

    if (amlt < 8.0) {
      dendamp = antp + damping * (amlt + 24.0 - phitp);
      switch0 = switchon(amlt, down_time - 24.0 - shift, sdel);
      geosync_trough = denmin * switch0 * (1.0 - switch1) + dendamp * (1.0 - switch0)
                       + dengrow * switch1 * (1.0 - switch2);
    } else {
      dendamp = antp + damping * (amlt - phitp);
      switch3 = switchon(amlt, down_time - shift, sdel);
      geosync_trough = denmin * switch3 + dengrow * switch1 * (1.0 - switch2)
                       + dendamp * switch2 * (1.0 - switch3);
    }

    // Scale density to L-shell of interest
    double trough_density = geosync_trough * pow(al, -4.5) / 2.0514092e-4;
    return trough_density;
  }

  /**
   * @brief This program calcuates the density profile for the interior
   * plasmasphere and for the emipirical modeling project. The Carpenter and
   * Anderson, 1992 plasmasphere is used only for a saturated plasmasphere. For a
   * more typical plasmaspheric density profile, we will use a profile obtained
   * from Gallagher et al. 1988
   * @param al (L-shell)
   * @param amlt (solar magnetic local time)
   * @param itime (time array)
   * @param am1 (output parameter)
   * @param b1 (output parameter)
   * @param x234 (output parameter)
   * @return double (inner plasmasphere density in cm^-3)
   */
  double ne_inner_ps(double al, double amlt, const std::array<int, 2>& itime, double& am1,
                     double& b1, double& x234) {
#if USE_STATIC_GCPM_CPP
    static int itime1_o = -99, itime2_o = -99;
    static double doy, doy_factor, rz12;
#else
    int itime1_o = -99, itime2_o = -99;
    double doy, doy_factor, rz12;
#endif

    const double a6 = -0.79, a7 = 5.208;
    double r = 1000.0 / RE + 1.0;

    if (itime[0] != itime1_o || itime[1] != itime2_o) {
      // call IRI to get the value of rz12 loaded
      IRIParams params = iri_sm(0.0, 0.0, r, itime);
      double rz12 = params.rz12;  // oarr[32]
      int iyear = itime[0] / 1000;
      doy = double(itime[0] - iyear * 1000);
      doy_factor = PI * (doy + 9.0) / 365.0;
      x234 = (0.15 * (cos(2.0 * doy_factor) - 0.5 * cos(4.0 * doy_factor))
              + (0.00127 * rz12 - 0.0635))
             * exp(-(al - 2.0) / 1.5);
      itime1_o = itime[0];
      itime2_o = itime[1];
    }

    //  calculate the inner plasmaspheric equatorial density
    // same calculation done in check_crossing
    // calculation used also in iri_ps_eq_bridge
    double inner_ps = std::pow(10.0, a6 * al + a7 + x234);

    am1 = a6;
    b1 = a7;

    return inner_ps;
  }

  /**
   * @brief c This routine checks to make sure that a8 is not beyond the point
   * where the inner plasmasphere and trough density models would cross were there
   * no plasmapause. If a8 is beyond that point, then the crossing point is
   * substituted for a8.
   *
   * @param a8
   * @param am1
   * @param b1
   * @param x234
   * @param amlt
   * @param akp
   * @return double
   */
  double check_crossing(double& a8, double am1, double b1, double x234, double amlt, double akp,
                        double geosync_trough) {
    // Determine where the inner plasmasphere plus plasmapause profile
    // intersects with the trough density profile.
    double stepl = 0.5;
    double zl = a8;
    double a = pow(10.0, am1 * zl + b1 + x234);
    double b = pp_profile(zl, amlt, akp, a8);
    double c = geosync_trough * pow(zl, -4.5) / 2.0514092e-4;
    double diff = a * b - c;
    int icount = 0;

    while (std::abs(stepl) > 0.05) {
      if ((diff < 0.0 && stepl > 0.0) || (diff > 0.0 && stepl < 0.0)) {
        stepl = -stepl / 2.0;
      }
      zl += stepl;
      // same calculation done in ne_inner_ps
      diff = std::pow(10, am1 * zl + b1 + x234) * pp_profile(zl, amlt, akp, a8)
             - geosync_trough * std::pow(zl, -4.5) / 2.0514092e-4;
      icount++;
      if (icount > 100) {
        std::cerr << "check_crossing reached loop-bound" << std::endl;
        return zl;
      }
    }
    return zl;
  }

  void iri_ps_eq_bridge(double al, double amlt, const std::array<int, 2>& itime, double& transh,
                        double& alpha, double& ano, double& am1, double& b1, double& x234,
                        double& psL) {
    double alatr = 0.0;
    double along = (amlt + 12.0) * AMLTRAD;
    along = along - (1.0 - sign(1.0, (12.0 - amlt))) * PI;

    IRIParams params;

    // get height and densiy of the f2 peak
    params = iri_sm(alatr, along, 200.0 / RE + 1.0, itime);
    double hmf2 = params.hmf2_km;
    double rf2 = hmf2 / RE + 1.0;

    // In an effort to reduce the cals to iri2007 the following is used
    // to approximate the point of maximum negative slope in the topside
    // ionosphere. This has been obtained from a linear fit to this location
    // (derived from the search algorithm above) as a function of returned
    // rz12 value from IRI2007. That analysis obtained this relationship:
    double rz12 = params.rz12;
    double ro_tmp = 1.05454 + 8.62678e-5 * rz12;
    double ro = std::max((rf2 + 0.01), ro_tmp);
    transh = (ro - 1.0) * RE;

    double ah1 = transh - 1.0;
    double ah2 = transh + 1.0;

    params = iri_sm(alatr, along, ro, itime);
    double dens = params.neiri;

    params = iri_sm(alatr, along, (ah1 / RE + 1.0), itime);
    double an1old = params.neiri;

    params = iri_sm(alatr, along, (ah2 / RE + 1.0), itime);
    double an2old = params.neiri;

    // calculate the initial equatorial power law transition function
    double alphao = -log10_safe(an1old / an2old) / log10_safe(ah1 / ah2);
    ano = dens / pow(transh, -alphao);

    // calculate where this initial power law intersects the plasmasphere profile
    double psh = 2000.0;
    for (int ii = 0; ii < 5; ++ii) {
      psh = pow(10.0, (am1 * (psh / RE + 1.0) + b1 + x234 + 6.0 - log10_safe(ano)) / (-alphao));
    }
    psL = psh / RE + 1.0;

    // Check to see whether psh is too large, meaning it was continously
    // increasing during the 5 iterations. If it is too large, then these two
    // curves did not intersect and we have to force a more steep topside
    // ionosphere.
    if (psh >= 0.5 * RE) {
      // Instead calculate the location along the equator where this power law
      // function matches the slope of the interior plasmaspheric density
      psL = 1.0 - alphao / am1 / log(10.0);
      psh = (psL - 1.0) * RE;
    }

    // new power law value, alpha, needs to match the ionosphere at the point
    // of maximum slope and the interior plasmaspheric density where the initial
    // power law slope matches the plasmasphere interior density slope
    double psden = pow(10.0, am1 * psL + b1 + x234 + 6.0);
    alpha = -log10_safe(dens / psden) / log10_safe(transh / psh);
    ano = dens / pow(transh, -alpha);
  }

  double ne_iri_cap(double r, double alatr, double amlt, const std::array<int, 2>& itime) {
    // Constants
    const double CHRRAD = 2.6179939e-1;
    const double AHCRIT = 350.0;
    const double OVERLAP = 50.0;
    const double POWERN = -2.8618;

    double aheight = (r - 1.0) * RE;
    double along = (amlt - 12.0) * CHRRAD;
    double ne_iri_cap_val;

    IRIParams params;

    if (aheight < (AHCRIT - OVERLAP)) {
      params = iri_sm(alatr, along, r, itime);
      ne_iri_cap_val = params.neiri;
    } else {
      params = iri_sm(alatr, along, (AHCRIT + RE) / RE, itime);
      double nb1 = params.neiri;
      double refn = std::log(nb1) + 16.764;

      ne_iri_cap_val = std::exp(POWERN * std::log(aheight) + refn) + 0.001;

      if (aheight <= (AHCRIT + OVERLAP)) {
        params = iri_sm(alatr, along, r, itime);
        double edensity_iri = params.neiri;

        double refhh = AHCRIT;
        double spred = -0.16;
        double widthh = OVERLAP;
        double refh2 = refhh - spred;
        double refh3 = refhh + spred;

        double switch2 = switchon(aheight, refh2, widthh);
        double switch3 = switchon(aheight, refh3, widthh);

        ne_iri_cap_val = edensity_iri * (1.0 - switch3) + ne_iri_cap_val * switch2;
      }
    }
    return ne_iri_cap_val;
  }

}  // namespace pecsim
