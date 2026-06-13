/**
 * @file gcpm/ne_iri.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-06
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <array>
#include <cmath>

#include "lupnt/environment/plasma/gcpm/iri_interface.h"

namespace pecsim {

  /**
   * @brief This subroutine is responsible for deteriming the total electron
   * density as a function of position within the ionosphere and plasmasphere.
   */
  double ne_iri_ps_trough(double r, double al, double alatr, double amlt, double akp,
                          const std::array<int, 2>& itime);

  /**
   * @brief This subroutine is responsible for deteriming the total electron
   * density as a function of position within the ionosphere and plasmasphere at
   * the magnetic equator.
   *
   */
  double ne_iri_ps_trough_eq(double al, double amlt, double akp, const std::array<int, 2>& itime);

  /**
   * @brief
  *  This subroutine models the plasmapause location and slope using a
  *  modified Lorentzian function as the basis.  It returns a factor
  *  that is used elsewhere to effect inclusion of the plasmapause
  *  transition from the inner plasmasphere to the trough.

  *  The factor returned varies from a value of 1 well inside a8
  *  to a value of 0 well outside a8.
  *  The rotation of the bulge with variation with Kpmax has been included
  *  A plasmapause profile is obtained from Carpenter & Anderson [1991] combined
  *  with Higel & Wu [1984] and Moldwin et al, [1994].
  *  The profile is included in a8.
  *  An average Kpmax dependence for the plasmapause slope is also included in
  *  a9.
   *
   * @param al
   * @param amlt  geomagnetic local time in hours
   * @param akp
   * @param a8
   * @return double
   */
  double pp_profile(double al, double amlt, double akp, double& a8);

  void bulge(double amlt, double akp, double& a8, double& a9, double& centroid);

  double ne_inner_ps(double al, double amlt, const std::array<int, 2>& itime, double& am1,
                     double& b1, double& x234);

  double ne_eq_trough(double al, double amlt, double akp, double& geosync_trough);

  double check_crossing(double& a8, double am1, double b1, double x234, double amlt, double akp,
                        double geosync_trough);

  /**
   * @brief Program runs the subroutine iri_sm to obtain IRI13 densities along an
   * L-shell.
   */
  void iri_ps_bridge(double r, double al, double alatr, double amlt,
                     const std::array<int, 2>& itime, double eq_iri_ps_trough, double& transh,
                     double& rf2, double& alpha, double& dno, double& co, double& switchh,
                     double& switchw, int& istat);

  /**
   * @brief Program determines the fit parameters for a power-law function
   * that is matched to the ionosphere at the point of maximum slope
   * and to the slope of the plasmaspheric interior density. This used
   * to be matched to the slope at the point of the maximum density
   * but that sometimes results in the power law function not falling
   * to the interior plasmaspheric densities. That won't work, so now
   * what is being done is to find the point of maximum slope, calculate
   * a power law function that matches the slope at that point, then
   * find where that function has the same slope as the interior plasmaspheric
   * density profile, then use that point and the point of maximum slope
   * to recalculate the power function. That will result in the power law
   * having a slope close to that of the topside ionosphere at the maximum slope
   * but still drop to interior plasmaspheric densities.
   *
   * First, then, find the location in the topside ionosphere where the slope is a
   * maximum. It starts looking at the F2 peak and stops when the slope starts to
   * decrease. It steps upward with an interval of delh.
   */
  void iri_ps_eq_bridge(double al, double amlt, const std::array<int, 2>& itime, double& transh,
                        double& alpha, double& ano, double& am1, double& b1, double& x234,
                        double& psL);

  /**
   * @brief This subroutine provides polar cap densities based on IRI and on the
   * GCPM polar cap model.
   *
   * @param r radial distance (RE)
   * @param alatr   geomagnetic latitude in radians
   * @param amlt    geomagnetic local time in hours
   * @param itime
   * @return double c The density bridge has the form:
   *    density = exp(-2.8618*log(r)+refn)
   *    where r is height in km
   *    this is obtained from the IRI density at 350km altitude and at the
   * latitude, L-shell, magnetic local time, and local time of the point provided
   * by the  empirical model user.  It is units of m-3, like the IRI total
   * electron density.
   *
   *    powern = -2.8618
   *    refalt = reference altitude where polar cap profile and IRI are made to
   * agree
   *
   *  The mathematical form of the density model comes from Persoon et al. and
   * Chandler et al, 1991. It approximates Alouette/ISIS, DE 1 RIMS, and DE 1 PWI
   * observations, while allowing for the density to rise and fall with the IRI
   * model through the refn parameter.
   */
  double ne_iri_cap(double r, double alatr, double amlt, const std::array<int, 2>& itime);

}  // namespace pecsim
