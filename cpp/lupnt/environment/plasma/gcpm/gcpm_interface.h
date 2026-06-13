/**
 * @file gcpm/gcpm_interface.h
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
#include <iostream>

#include "lupnt/environment/plasma/core/math_utils.h"
#include "lupnt/environment/plasma/core/user_filepath.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"
#include "lupnt/environment/plasma/gcpm/conversions.h"
#include "lupnt/environment/plasma/gcpm/iri_interface.h"
#include "lupnt/environment/plasma/gcpm/ne_iri.h"

// Placeholder for external function implementation
namespace pecsim {

  /**
   * @brief c		Global Core Plasma Model ***Version 2.3***
   * Original Version 1.0, January 1, 2000
   *
   * modified by dlg 8/3/2007 to include the polar cap model
   * modified by dlg 1/6/2009 to include seasonal and solar cycle variations in
   * inner plasmasphere from C&A 1992 corrected by dlg 6/3/2009 to fix
   * ionosphere-plasmasphere bridge code that was trying to make bridge below the
   * F2 peak.
   *
   * @param datatime  date and time (in UTC)
   * @param r      geocentric radial distance in Re
   * @param amlt   solar magnetic local time in hours
   * @param alatr  solar magnetic latitude in radians
   * @param akp    planetary Kp index (-1: use the Kp index from data)
   * @return std::vector<double>
   *      [0] = total electron density in 1/cm^3
   * 	    [1] = total hydrogen density in 1/cm^3
   *      [2] = total helium density in 1/cm^3
   *      [3] = total oxygen density in 1/cm^3]
   *      [4] = f107 index
   *      [5] = rz12 index
   *      [6] = hmf2_km (F2 peak height in km)
   *
   */
  std::vector<double> gcpm_v24(DateTime datetime, double r, double amlt, double alatr,
                               double akp = -1);

  /**
   * @brief Call the GCPM Fortran function
   *
   * @param datetime  The date and time in UTC
   * @param r_RE    The geocentric distance in Earth radii
   * @param amlt    The solar magnetic local time in hours
   * @param alatr   The solar magnetic latitude in radians
   * @param akp     The KP index (-1: use the KP index from data)
   * @return std::vector<double>   The output density array
   * (4 elements, total electron density, total hydrogen density, helium density,
   * oxygen density) in cm^-3
   */
  std::vector<double> gcpm_v24_fortran(DateTime datetime, double r_RE, double amlt, double alatr,
                                       double akp = -1);
}  // namespace pecsim
