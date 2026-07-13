/**
 * @file gcpm_interface.cpp
 * @author Keidai Iiyama
 * @brief  This file contains the implementation of the GCPM model interface
 * @version 0.1
 * @date 2025-02-06
 *
 * @copyright Copyright (c) 2025
 *
 */
#include "lupnt/environment/plasma/gcpm/gcpm_interface.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/core/math_utils.h"
#include "lupnt/environment/plasma/core/user_filepath.h"
#include "lupnt/environment/plasma/env/kp.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {

  extern "C" {
  void gcpm_v24_(int* itime, float* r, float* amlt, float* alatr, float* akp, float* outn);
  }

  std::vector<double> gcpm_v24(DateTime datetime, double r, double amlt, double alatr, double akp) {
    // first move to the data directory
    std::string curr_dir = std::filesystem::current_path();
    std::filesystem::current_path(get_iri_data_dir());

    // get kp index
    if (akp < 0) {
      akp = get_kp_index(datetime);
    }

    // Convert to itime
    std::array<int, 2> itime;
    itime[0] = datetime.year * 1000 + datetime.doy;  // Year and day of the year
    int ihr = datetime.hour, imin = datetime.min, isec = datetime.sec;
    itime[1] = int((ihr * 3600 + imin * 60 + isec) * 1000);  // Milliseconds of the day

    // Constants
    const double re = RE;  // Earth radius in km
    double clat, al, aheight;

    // the half width in L-shell of over which the transition
    // takes place between the trough and polar cap models.
    // The L-shell at which this transition is centered is
    // obtained from PN, which hold empirical locations for
    // the polarward edge of the auroral zone for several
    // values of Kp and magnetic local time.
    double altrans = 2.0;

    // Get f107 value from IRI
    IRIParams params = iri_sm(0.0, 0.0, 1.01, itime);
    double f107 = params.f107;
    double rz12 = params.rz12;
    double hmf2_km = params.hmf2_km;  // F2 peak height in km

// Persistent variables
#if USE_STATIC_GCPM_CPP
    static double oldmlt = -1.0, oldkp = -1.0;
    static double tranlow = 0.0, tranhigh = 0.0;
    static double alcrit;
#else
    double oldmlt = -1.0, oldkp = -1.0;
    double tranlow = 0.0, tranhigh = 0.0;
    double alcrit;
#endif

    // Will execute this section if we need new a new value for the
    // invarient latitude of the polarward edge of the auroral zone.
    // This location is determined as a function of MLT and Kp from
    // the array PN.
    oldmlt = amlt;
    oldkp = akp;
    double bmlt = amlt * 3.0 + 1.0;

    int imlt = std::floor(bmlt);
    double diffmlt = bmlt - static_cast<double>(imlt);
    if (imlt > 72) imlt = 1;
    int jmlt = imlt + 1;
    if (jmlt > 72) jmlt = 1;

    int ikp = static_cast<int>(akp + 1.0);
    double diffkp = akp - std::floor(akp);
    if (ikp > 10) ikp = 10;
    int jkp = ikp + 1;
    if (jkp > 10) jkp = 10;

    // Fortran indexing to C++ indexing
    imlt = imlt - 1;
    jmlt = jmlt - 1;
    ikp = ikp - 1;
    jkp = jkp - 1;

    double pn1 = (PN(jmlt, ikp) - PN(imlt, ikp)) * diffmlt + PN(imlt, ikp);
    double pn2 = (PN(jmlt, jkp) - PN(imlt, jkp)) * diffmlt + PN(imlt, jkp);
    double ps1 = (PS(jmlt, ikp) - PS(imlt, ikp)) * diffmlt + PS(imlt, ikp);
    double ps2 = (PS(jmlt, jkp) - PS(imlt, jkp)) * diffmlt + PS(imlt, jkp);

    double alatcrits = (ps2 - ps1) * diffkp + ps1;
    double alatcritn = (pn2 - pn1) * diffkp + pn1;

    alcrit = 1.0 / pow(cos(alatcritn * DEG2RAD), 2);
    double aurora_mlat = alatcritn;
    tranlow = alcrit - altrans;
    tranhigh = alcrit + altrans;

    // We need to obtain the L-shell of the location given, while limiting
    // the maximum L-shell that will be used.  Higher latitudes and L-shells
    // technically extending to infinity will be described by the maximum
    // acceptable L-shell value, without causing problems.
    clat = pow(cos(alatr), 2);
    if (clat < 1.0e-5) clat = 1.0e-5;
    al = r / clat;
    aheight = (r - 1.0) * re;

    // if (r > 1.0485) {
    //     std::cout << "altitude high" << std::endl;
    // }

    // The model contains elements for the plasmasphere and plasmapause, the
    // trough, and the polar cap.  The IRI is used for the ionosphere.
    // The polar cap - trough boundary is chosen to be at L=alcrit.  The function
    // subroutine switchon is used to transition between these regions.  This
    // function is also used to transition between the IRI and the other exterior
    // models.
    double edensity = 0.0;
    std::string call_type = "";
    if (al < tranlow) {
      // execute this section if we are equatorward of the polar cap
      edensity = ne_iri_ps_trough(r, al, alatr, amlt, akp, itime);
      call_type = "low";
    } else if (al <= tranhigh) {
      // execute this section if we are in the transition region between the
      // trough and polar cap regions
      double ps_edensity = ne_iri_ps_trough(r, al, alatr, amlt, akp, itime);
      double cap_edensity = ne_iri_cap(r, alatr, amlt, itime);
      double switch_val = switchon(al, alcrit, altrans);
      edensity = ps_edensity * (1.0 - switch_val) + cap_edensity * switch_val;
      call_type = "mid";
    } else {
      // execute this section if we are polarward of the trough
      edensity = ne_iri_cap(r, alatr, amlt, itime);
      call_type = "polar";
    }

    // Convert to particles per cm^3
    double den = edensity / 1.0e6;

    // Compute He+ to H+ density ratio in the plasmasphere
    double aHeH = pow(10.0, (-1.541 - 0.176 * r + 8.557e-3 * f107 - 1.458e-5 * f107 * f107));

    // Helium concentration drops dramatically with transition to high latitudes
    // and open field lines
    aHeH *= (1.0 - switchon(al, (tranlow + tranhigh) / 2.0, altrans));

    // Compute relative O+ density
    double alphaO = 0.995 / pow((1.0 + pow((aheight - 350.0) / 531.25, 2)), 3) + 0.005;

    // Compute relative He+ concentration in the plasmasphere
    double alphaHe = (aHeH != 0.0) ? (1.0 - alphaO) / (1.0 + 1.0 / aHeH) : 0.0;

    // Compute densities of H+, He+, and O+
    double heden = alphaHe * den;        // Helium density [2]
    double oxden = alphaO * den;         // Oxygen density  [3]
    double hyden = den - heden - oxden;  // Hydrogen density [1]

    double call_type_idx = 0;
    if (call_type == "low") {
      call_type_idx = 0;  // Low latitude
    } else if (call_type == "mid") {
      call_type_idx = 1;  // Mid latitude
    } else if (call_type == "polar") {
      call_type_idx = 2;  // Polar latitude
    }

    std::vector<double> outn = {den, hyden, heden, oxden, f107, rz12, hmf2_km, call_type_idx};

    // move back to the original directory
    std::filesystem::current_path(curr_dir);

    return outn;
  }

  std::vector<double> gcpm_v24_fortran(DateTime datetime, double r_RE, double amlt, double alatr,
                                       double akp) {
    // Convert the input arguments to float
    // first move to the data directory
    std::string curr_dir = std::filesystem::current_path();
    std::filesystem::current_path(get_iri_data_dir());

    // get kp index
    if (akp < 0) {
      akp = get_kp_index(datetime);
    }

    //
    int itime[2];
    itime[0] = datetime.year * 1000 + datetime.doy;  // Year and day of the year
    int ihr = datetime.hour, imin = datetime.min, isec = datetime.sec;
    itime[1] = (ihr * 3600 + imin * 60 + isec) * 1000;  // Milliseconds of the day
    float r_f = static_cast<float>(r_RE);               // Geo-centric distance in Earth radii
    float amlt_f = static_cast<float>(amlt);            // Magnetic local time in hours
    float alatr_f = static_cast<float>(alatr);          // Magnetic latitude in radians
    float akp_f = static_cast<float>(akp);              // KP index

    float outn[4];

    // Call the Fortran function
    gcpm_v24_(itime, &r_f, &amlt_f, &alatr_f, &akp_f, outn);

    // Same future-epoch robustness as iri_sm(): the shared IRI solar-index data
    // (apf107.dat / ig_rz.dat) ends ~2024 with only a short projection, so for
    // later epochs the Fortran GCPM returns NaN densities. Step the requested
    // year back (itime[0] = year*1000 + doy) until it has solar data, so a
    // future-epoch call falls back to the most recent activity instead of NaN.
    for (int back = 1; back <= 30 && !(outn[0] >= 0.0f); ++back) {
      int itime_fb[2] = {itime[0] - back * 1000, itime[1]};
      gcpm_v24_(itime_fb, &r_f, &amlt_f, &alatr_f, &akp_f, outn);
    }

    std::vector<double> outvec;
    for (int i = 0; i < 4; ++i) {
      outvec.push_back(static_cast<double>(outn[i]));
    }

    // move back to the original directory
    std::filesystem::current_path(curr_dir);

    return outvec;
  }

}  // namespace pecsim
