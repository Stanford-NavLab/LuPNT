/**
 * @file iri_model.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-07
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "lupnt/environment/plasma/gcpm/iri_interface.h"

#include <omp.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/core/user_filepath.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"
#include "lupnt/environment/plasma/gcpm/conversions.h"

namespace pecsim {

  extern "C" {
  // 2007 version
  void initialize_();
  void readapf107_2007_(void);
  void iri_sub_(int* jf, int* jmag, float* blatd, float* blongd, int* yyyy, int* ddd, float* dhour,
                float* aheight1, float* aheight2, float* delh, float* outf, float* oarr);

  // 2020 version
  void read_ig_rz_2020_(void);
  void readapf107_2020_(void);
  void iri_sub_2020_(int* jf, int* jmag, float* blatd, float* blongd, int* yyyy, int* ddd,
                     float* dhour, float* aheight1, float* aheight2, float* delh, float* outf,
                     float* oarr);
  }

  // IRI Model selector ---------------------------------
  static IRIModel iri_model = IRIModel::NONE;
  static IRI2007Option iri_option_2007 = IRI2007Option();
  static IRI2020Option iri_option_2020 = IRI2020Option();

  void print_outf(const float* outf, int rows, int cols) {
    std::cout << "outf" << std::endl;
    for (int i = 0; i < rows; i++) {
      std::cout << "Row " << i << ": ";
      for (int j = 0; j < cols; j++) {
        std::cout << outf[i + j * rows] << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  void print_oarr(const float* oarr, int rows, int cols) {
    std::cout << "oarr" << std::endl;
    int k = 0;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        std::cout << k + 1 << " : " << oarr[k] << ",  ";
        k++;
      }
      std::cout << std::endl;
    }
  }

  void set_iri_model(IRIModel model) {
    iri_model = model;
    std::filesystem::current_path(get_iri_data_dir());

#pragma omp critical(file_io)
    {
      if (model == IRIModel::IRI_2007) {
        readapf107_2007_();
      } else if (model == IRIModel::IRI_2020) {
        readapf107_2020_();
        read_ig_rz_2020_();
      }
    }
  }

  IRIModel get_iri_model() {
    if (iri_model == IRIModel::NONE) {
      std::runtime_error("IRI model not set.");
    }
    return iri_model;
  }

  std::string get_iri_model_str() {
    if (iri_model == IRIModel::IRI_2007) {
      return "IRI 2007";
    } else if (iri_model == IRIModel::IRI_2020) {
      return "IRI 2020";
    } else {
      return "Unknown IRI Model";
    }
  }

  std::string get_iri_data_dir() { return get_base_path() + "/data/iri"; }

  void set_iri2007_option(const IRI2007Option& option) {
    iri_option_2007.R12 = option.R12;
    iri_option_2007.update_jf_2007();  // Update the Fortran options
  }

  void set_iri2020_option(const IRI2020Option& option) {
    iri_option_2020.compute_teti = option.compute_teti;
    iri_option_2020.compute_ni = option.compute_ni;
    iri_option_2020.output_messages = option.output_messages;
    iri_option_2020.output_to_text = option.output_to_text;
    iri_option_2020.plasma_model = option.plasma_model;
    iri_option_2020.without_plasmapause = option.without_plasmapause;
    iri_option_2020.R12 = option.R12;
    iri_option_2020.update_jf_2020();  // Update the Fortran options
  }

  void IRI2007Option::update_jf_2007() {
    int jf[30] = {true,  true,  true,  true, false, false, true, true,  true,  true,
                  true,  false, true,  true, true,  true,  true, true,  true,  true,
                  false, true,  false, true, true,  true,  true, false, false, false};

    // set options based on the IRI2007Option
    if (R12 > 0 && R12 <= 200) {
      jf[25] = 0; /* storm model off */
      jf[16] = 0; /* user input R12 */
      jf[24] = 0; /* user input F10.7 */
      jf[26] = 0; /* user input IG12 */
    } else if (R12 == -1) {
      jf[25] = 1; /* storm model on */
      jf[16] = 1; /* historical or projected R12 */
      jf[24] = 1; /* historical or projected F10.7 */
      jf[26] = 1; /* historical of projected IG12 */
    } else if (R12 == -2) {
      jf[25] = 0; /* storm model off */
      jf[16] = 1; /* historical or projected R12 */
      jf[24] = 1; /* historical or projected F10.7 */
      jf[26] = 1; /* historical of projected IG12 */
    } else {
      throw std::runtime_error("Invalid value for R12");
    }

    for (int i = 0; i < 30; i++) {
      jf2007[i] = jf[i];  // initialize all to true
    }
  }

  void IRI2020Option::update_jf_2020() {
    std::vector<int> false_lists = {4, 5, 6, 23, 30, 33, 35, 39, 40, 47};  // false index in fortran

    int jf[50];

    for (int i = 0; i < 50; i++) {
      jf[i] = true;  // initialize all to true
    }
    for (int i = 0; i < false_lists.size(); i++) {
      jf[false_lists[i] - 1] = false;
    }

    // set options based on the IRI2020Option
    jf[1 - 1] = compute_teti;          // compute TETI
    jf[2 - 1] = compute_ni;            // compute NI
    jf[12 - 1] = 1 - output_to_text;   // output to text file (1: yes, 0: no)
    jf[34 - 1] = output_messages;      // turn off messages
    jf[49 - 1] = plasma_model;         // plasma model (1: ozhogin, 0: gallager)
    jf[50 - 1] = without_plasmapause;  // without plasmapause

    // set options based on the IRI2007Option
    if (R12 > 0 && R12 <= 200) {
      jf[16] = 0; /* user input R12 */
      jf[24] = 0; /* user input F10.7 */
      jf[26] = 0; /* user input IG12 */
      jf[31] = 0; /* user input F10.7_81 */
      jf[25] = 0; /* default setting for foF2 storm model is off */
    } else if (R12 == -1.0) {
      jf[16] = 1; /* historical or projected R12 */
      jf[24] = 1; /* historical or projected F10.7 */
      jf[26] = 1; /* historical or projected IG12 */
      jf[31] = 1; /* historical or projected F10.7_81 */
      jf[25] = 1; /* default setting for foF2 storm model is on */
    } else {
      throw std::runtime_error("Invalid value for R12");
    }

    for (int i = 0; i < 50; i++) {
      jf2020[i] = jf[i];  // initialize all to true
    }
  }

  std::vector<double> iri_2020(DateTime datetime, double r, double amlt, double alatr, double akp,
                               IRI2020Option op) {
    // set version to IRI 2020
    set_iri_model(IRIModel::IRI_2020);
    set_iri2020_option(op);  // Set the IRI options

    // Convert to itime
    std::array<int, 2> itime = datetime_to_itime(datetime);
    double along = lt_to_long(amlt);  // Convert magnetic local time to longitude

    // Get values value from IRI
    IRIParams params = iri_sm(alatr, along, r, itime);

    // Prepare output vector
    std::vector<double> outn(4);
    outn[0] = params.neiri * 1e-6;   // Total electron density
    outn[1] = params.nhoiri * 1e-6;  // H+ density
    outn[2] = params.nheiri * 1e-6;  // He+ density
    outn[3] = params.noiri * 1e-6;   // O+ density

    return outn;
  }

  std::vector<double> iri_2007(DateTime datetime, double r, double amlt, double alatr, double akp,
                               IRI2007Option op) {
    // set version to IRI 2007
    set_iri_model(IRIModel::IRI_2007);
    set_iri2007_option(op);  // Set the IRI options

    // Convert to itime
    std::array<int, 2> itime = datetime_to_itime(datetime);
    double along = lt_to_long(amlt);  // Convert magnetic local time to longitude

    // Get values value from IRI
    IRIParams params = iri_sm(alatr, along, r, itime);

    // Prepare output vector
    std::vector<double> outn(4);
    outn[0] = params.neiri * 1e-6;   // Total electron density
    outn[1] = params.nhoiri * 1e-6;  // H+ density
    outn[2] = params.nheiri * 1e-6;  // He+ density
    outn[3] = params.noiri * 1e-6;   // O+ density

    return outn;
  }

  IRIParams iri_sm(double alatr, double along, double r, const std::array<int, 2>& itime) {
    int yyyy = itime[0] / 1000;
    int ddd = itime[0] - yyyy * 1000;
    double dhour = static_cast<double>(itime[1]) / 3600000.0 + 25.0;  // convert to UT

    double aheight = (r - 1.0) * RE;

    IRIParams params = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (aheight > 3000.0) {
      return params;
    }

    std::array<double, 3> pos_sm, pos_geo;
    pol_to_cart(alatr, along, r, pos_sm);
    sm_to_geo(itime, pos_sm, pos_geo);

    double blatr, blong, rtemp;
    cart_to_pol(pos_geo, blatr, blong, rtemp);
    double blatd = blatr * RAD2DEG;
    double blongd = blong * RAD2DEG;

    int jmag = 0;  // use geographic coordinates

    params = iri_sub(jmag, blatd, blongd, yyyy, -ddd, dhour, aheight, aheight, DELH);

    return params;
  }

  // This subroutine is used to call the IRI model from GCPM
  IRIParams iri_sub(int jmag, double blatd, double blongd, int yyyy, int mddd, double dhour,
                    double aheight1, double aheight2, double delh) {
    // Placeholder for IRI model subroutine

    double rz12, f107, hmf2_km, ne, nh, nhe, no;

    IRIModel model = get_iri_model();

    // set parameters for IRI model
    float blatdf = static_cast<float>(blatd);
    float blongdf = static_cast<float>(blongd);
    float dhourf = static_cast<float>(dhour);
    float aheight1f = static_cast<float>(aheight1);
    float aheight2f = static_cast<float>(aheight2);
    float delhf = static_cast<float>(delh);

    double IG12, R12;

    // IRI 2007 model -------------------------------------------------------
    if (model == IRIModel::IRI_2007) {
      float outf[20 * 100] = {0};
      float oarr[50] = {0};

      // Calculate IG12 and F10.7 from the user input R12
      R12 = iri_option_2007.R12;
      if (R12 > 0) {
        f107 = 63.75 + R12 * (0.728 + R12 * 0.00089);
        IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
        *(oarr + 32) = static_cast<float>(R12);
        *(oarr + 38) = static_cast<float>(IG12);
        *(oarr + 40) = static_cast<float>(f107);
      }

      iri_sub_(iri_option_2007.jf2007, &jmag, &blatdf, &blongdf, &yyyy, &mddd, &dhourf, &aheight1f,
               &aheight2f, &delhf, outf, oarr);

      // print_outf(outf, 20, 25);
      // print_oarr(oarr, 5, 10);

      outf[0] = static_cast<float>(std::max(0.0, static_cast<double>(outf[0])));

      rz12 = static_cast<double>(oarr[32]);
      f107 = static_cast<double>(oarr[40]);
      hmf2_km = static_cast<double>(oarr[1]);
      ne = static_cast<double>(outf[0]);
      nh = static_cast<double>(outf[5]);
      nhe = static_cast<double>(outf[6]);
      no = static_cast<double>(outf[4]);
    }
    // IRI 2020 model -------------------------------------------------------
    else if (model == IRIModel::IRI_2020) {
      float outf[20 * 1000] = {0.0f};
      float oarr[100] = {0.0f};

      for (int i = 0; i < 20 * 1000; i++) {
        outf[i] = 0.0f;  // initialize all to 0.0
      }
      for (int i = 0; i < 100; i++) {
        oarr[i] = -1.0f;  // initialize all to 0.0
      }

      // Calculate IG12 and F10.7 from the user input R12
      R12 = iri_option_2020.R12;
      if (R12 > 0.0) {
        f107 = 63.75 + R12 * (0.728 + R12 * 0.00089);
        IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
        *(oarr + 32) = static_cast<float>(R12);
        *(oarr + 38) = static_cast<float>(IG12);
        *(oarr + 40) = static_cast<float>(f107);
        *(oarr + 45) = static_cast<float>(f107); /* F10.7_81 is set to F10.7 */
      }

      // Call the IRI2020 subroutine
      iri_sub_2020_(iri_option_2020.jf2020, &jmag, &blatdf, &blongdf, &yyyy, &mddd, &dhourf,
                    &aheight1f, &aheight2f, &delhf, outf, oarr);

      outf[0] = static_cast<float>(std::max(0.0, static_cast<double>(outf[0])));

      rz12 = static_cast<double>(oarr[32]);
      f107 = static_cast<double>(oarr[40]);
      hmf2_km = static_cast<double>(oarr[1]);
      ne = static_cast<double>(outf[0]);
      nh = static_cast<double>(outf[5]);
      nhe = static_cast<double>(outf[6]);
      no = static_cast<double>(outf[4]);
    } else {
      std::runtime_error("IRI model not set.");
    }

    // These values can be used further as needed
    IRIParams params = {rz12, f107, hmf2_km, ne, nh, nhe, no};

    return params;
  }  // namespace pecsim

}  // namespace pecsim
