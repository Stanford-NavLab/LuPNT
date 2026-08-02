/**
 * @file gcpm/iri_interface.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-07
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {
  enum IRIModel { NONE, IRI_2007, IRI_2020 };

  struct IRIParams {
    double rz12;     // oarr[32] (oarr[33] in fortran)
    double f107;     // oarr[40] (oarr[41] in fortran)
    double hmf2_km;  // oarr[1] (oarr[2] in fortran)  F2 peak height in km
    double neiri;    // outf[0][0]
    double nhoiri;   // outf[5][0]
    double nheiri;   // outf[6][0]
    double noiri;    // outf[4][0]
  };

  class IRI2007Option {
  public:
    double R12 = -1.0;  // R12 index, >0: use value as R12,
                        // -1: use historical or projected R12 (with storm model)
                        // -2: no storm model, use historical or projected R12

    // GCPM consumes only the IRI electron density, hmF2, and F10.7/R12, never
    // the ion temperatures or ion composition. These default to false so the
    // GCPM density path (and compute_ne) skips those expensive IRI sub-models;
    // the standalone iri_2007() accessor forces them back on where its returned
    // ion densities require them.
    bool compute_teti = false;  // jf(2): electron/ion temperatures
    bool compute_ni = false;    // jf(3): ion composition (H+/He+/O+)

    int jf2007[30];  // Fortran IRI 2007 options

    IRI2007Option() { update_jf_2007(); }

    void update_jf_2007();
  };

  class IRI2020Option {
  public:
    bool compute_teti = true;      // 1
    bool compute_ni = true;        // 2
    bool output_messages = false;  // 34
    bool output_to_text = false;   // 12 (1: yes, 0: no)
    int plasma_model = 1;          // 49 (1: ozhogin  0: gallager)
    int without_plasmapause = 1;   // 50
    double R12 = -1.0;             // R12 index, not used in IRI2020

    int jf2020[50];  // Fortran IRI 2020 options

    IRI2020Option() { update_jf_2020(); }

    void update_jf_2020();
  };

  // IRI model selection
  void set_iri_model(IRIModel model);
  IRIModel get_iri_model();
  std::string get_iri_model_str();
  std::string get_iri_data_dir();
  void set_iri2007_option(const IRI2007Option& option);
  void set_iri2020_option(const IRI2020Option& option);

  // Call the IRI model from the same input as gcpm_v24
  std::vector<double> iri_2020(DateTime datetime, double r, double amlt, double alatr, double akp,
                               IRI2020Option op = IRI2020Option());
  std::vector<double> iri_2007(DateTime datetime, double r, double amlt, double alatr, double akp,
                               IRI2007Option op = IRI2007Option());

  // This routine is used to call the IRI model from GCPM
  IRIParams iri_sm(double alatr, double along, double r, const std::array<int, 2>& itime);

  // Steps the year of `itime` back to the most recent one the bundled IRI solar-index data
  // (apf107.dat / ig_rz.dat) covers, so downstream GCPM/IRI density paths return finite
  // values instead of the no-data NaN sentinel for future epochs. Self-adjusting (no
  // hard-coded cutoff year); returns `itime` unchanged when it is already covered. Assumes
  // the IRI data directory is the current working path (as gcpm_v24 sets it).
  std::array<int, 2> iri_valid_solar_itime(const std::array<int, 2>& itime);

  // This subroutine is used to call the IRI model from GCPM
  IRIParams iri_sub(int jmag, double blatd, double blongd, int yyyy, int ddd, double dhour,
                    double aheight1, double aheight2, double delh);

}  // namespace pecsim
