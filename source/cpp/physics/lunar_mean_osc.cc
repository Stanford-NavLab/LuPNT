/**
 * @file LunarMeanOsc.cpp
 * @author Stanford NAV LAB
 * @brief  Lunar Mean and Osculating Elements
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <array>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  Vec6 Mean2Osculatingulating(Vec6 meanCoeVec) {
    //  double t = 0.0;  // Time [s]
    //  double nM = 2.66e-6;
    Vec6 meanDoeVec = Classical2Delaunay(meanCoeVec, GM_MOON);

    double lpp = (double)meanDoeVec(0);
    double gpp = (double)meanDoeVec(1);
    double hpp = (double)meanDoeVec(2);
    double Lpp = (double)meanDoeVec(3);
    double Gpp = (double)meanDoeVec[4];
    double Hpp = (double)meanDoeVec[5];

    std::array<double, 6> sp2 = ComputeSecondOrderShortPeriod(meanCoeVec, meanDoeVec);
    std::array<double, 6> mp1 = ComputeFirstOrderMediumPeriod(meanCoeVec, meanDoeVec);
    std::array<double, 6> mp2 = ComputeSecondOrderMediumPeriod(meanCoeVec, meanDoeVec);

    double l = lpp + sp2[0] + mp1[0] + mp2[0];  // l = lpp + lpp_sp2 + lpp_mp1 + lpp_mp2;
    double g = gpp + sp2[1] + mp1[1] + mp2[1];  // g = gpp + gpp_sp2 + gpp_mp1 + gpp_mp2;
    double h = hpp + sp2[2] + mp1[2] + mp2[2];  // h = hpp + hpp_sp2 + hpp_mp1 + hpp_mp2;
    double L = Lpp + sp2[3] + mp1[3] + mp2[3];  // L = Lpp + Lpp_sp2 + Lpp_mp1 + Lpp_mp2;
    double G = Gpp + sp2[4] + mp1[4] + mp2[4];  // G = Gpp + Gpp_sp2 + Gpp_mp1 + Gpp_mp2;
    double H = Hpp + sp2[5] + mp1[5] + mp2[5];  // H = Hpp + Hpp_sp2 + Hpp_mp1 + Hpp_mp2;

    Vec6 oscDoeVec;
    oscDoeVec << l, g, h, L, G, H;

    Vec6 oscCoeVec = Delaunay2Classical(oscDoeVec, GM_MOON);
    return oscCoeVec;
  }

  Vec6 Osculating2Mean(Vec6 oscCoeVec) {
    //  double t = 0.0;
    //  double nM = 2.66e-6;
    Vec6 oscDoeVec = Classical2Delaunay(oscCoeVec, GM_MOON);

    double l = (double)oscDoeVec(0);
    double g = (double)oscDoeVec(1);
    double h = (double)oscDoeVec(2);
    double L = (double)oscDoeVec(3);
    double G = (double)oscDoeVec[4];
    double H = (double)oscDoeVec[5];

    std::array<double, 6> sp2 = ComputeSecondOrderShortPeriod(oscDoeVec, oscDoeVec);
    std::array<double, 6> mp1 = ComputeFirstOrderMediumPeriod(oscDoeVec, oscDoeVec);
    std::array<double, 6> mp2 = ComputeSecondOrderMediumPeriod(oscDoeVec, oscDoeVec);
    std::array<double, 6> mp1c = ComputeCorrectionMediumPeriod(oscCoeVec, oscDoeVec);

    double lpp
        = l - sp2[0] - mp1[0] - mp1c[0] - mp2[0];  // lpp = l - l_sp2 - l_mp1 - l_mp1c - l_mp2;
    double gpp
        = g - sp2[1] - mp1[1] - mp1c[1] - mp2[1];  // gpp = g - g_sp2 - g_mp1 - g_mp1c - g_mp2;
    double hpp
        = h - sp2[2] - mp1[2] - mp1c[2] - mp2[2];  // hpp = h - h_sp2 - h_mp1 - h_mp1c - h_mp2;
    double Lpp
        = L - sp2[3] - mp1[3] - mp1c[3] - mp2[3];  // Lpp = L - L_sp2 - L_mp1 - L_mp1c - L_mp2;
    double Gpp
        = G - sp2[4] - mp1[4] - mp1c[4] - mp2[4];  // Gpp = G - G_sp2 - G_mp1 - G_mp1c - G_mp2;
    double Hpp
        = H - sp2[5] - mp1[5] - mp1c[5] - mp2[5];  // Hpp = H - H_sp2 - H_mp1 - H_mp1c - H_mp2;

    Vec6 meanDoeVec;
    meanDoeVec << lpp, gpp, hpp, Lpp, Gpp, Hpp;

    Vec6 meanCoeVec = Delaunay2Classical(meanDoeVec, GM_MOON);
    return meanCoeVec;
  }

  std::array<double, 6> ComputeSecondOrderShortPeriod(Vec6 &coe, Vec6 &doe) {
    double n3 = 2.66e-6;
    //  double nM = 2.66e-6;
    //  double J2 = 2.03e-4;
    double k = 0.98785;

    double a = (double)coe(0);
    double e = (double)coe(1);
    double i = (double)coe(2);
    //  double O = (double)coe(3);
    //  double w = (double)coe(4);
    //  double M = (double)coe(5);
    double f = (double)Mean2TrueAnomaly(coe(5), coe(1));

    double l = (double)doe(0);
    double g = (double)doe(1);
    double h = (double)doe(2);
    //  double L = (double)doe(3);
    //  double G = (double)doe(4);
    //  double H = (double)doe(5);

    // The second-order terms of the short-period variations are

    // The coefficients in the generating function S∗ 2 are given by (Giacaglia et
    // al. 1970):
    double B_2_1 = e * sin(f) + (f - l);
    double B_2_2 = 3 * e * sin(f + 2 * g) + 3 * sin(2 * f + 2 * g) + e * sin(3 * f + 2 * g);
    double B_22_1 = 2 * (f - l) * cos(2 * h) + e * sin(f - 2 * h) + e * sin(f + 2 * h);
    double B_22_2 = 3 * e * sin(f + 2 * g - 2 * h) + 3 * sin(2 * f + 2 * g - 2 * h)
                    + e * sin(3 * f + 2 * g - 2 * h);
    double B_22_3 = 3 * e * sin(f + 2 * g + 2 * h) + 3 * sin(2 * f + 2 * g + 2 * h)
                    + e * sin(3 * f + 2 * g + 2 * h);
    double B_M_1
        = 9 * e * (4 + pow(e, 2)) * sin(E) - 9 * pow(e, 2) * sin(2 * E) + pow(e, 3) * sin(3 * E);
    double B_M_2 = 2 * cos(2 * h)
                       * (-9 * e * (4 + pow(e, 2)) * sin(E) + 9 * pow(e, 2) * sin(2 * E)
                          - pow(e, 3) * sin(3 * E))
                   + 2 * cos(2 * g)
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E))
                   + 2 * sin(2 * g)
                         * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                            + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                            - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E));
    double B_M_3 = 2 * sin(2 * g - 2 * h)
                       * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                          + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                          - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                   + 2 * cos(2 * g - 2 * h)
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E));
    double B_M_4 = 2 * sin(2 * g + 2 * h)
                       * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                          + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                          - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                   + 2 * cos(2 * g + 2 * h)
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E));

    // B_M_e
    double B_M_e_0 = 6 * sin(E) * (6 + 5 * pow(e, 2) - 6 * e * cos(E) + pow(e, 2) * cos(2 * E));
    double B_M_e_1
        = -12 * sin(E) * (6 + 5 * pow(e, 2) - 6 * e * cos(E) + pow(e, 2) * cos(2 * E)) * cos(2 * h)
          + 4 * sin(E) * ((3 * pow(e, 2) - 2) * cos(2 * E) - 21 * pow(e, 2) + 6 * e * cos(E) - 16)
                * cos(2 * g)
          + (4 / sqrt(1 - pow(e, 2)))
                * (15 * (2 * pow(e, 2) - 1) * cos(E) + 3 * e * (1 - 3 * pow(e, 2)) * cos(2 * E)
                   + (2 * pow(e, 2) - 1) * cos(3 * E))
                * sin(2 * g);
    double B_M_e_2
        = 4 * sin(E) * ((3 * pow(e, 2) - 2) * cos(2 * E) - 21 * pow(e, 2) + 6 * e * cos(E) - 16)
              * cos(2 * g - 2 * h)
          + (4 / sqrt(1 - pow(e, 2)))
                * (15 * (2 * pow(e, 2) - 1) * cos(E) + 3 * e * (1 - 3 * pow(e, 2)) * cos(2 * E))
                * sin(2 * g - 2 * h);
    double B_M_e_3
        = 4 * sin(E) * ((3 * pow(e, 2) - 2) * cos(2 * E) - 21 * pow(e, 2) + 6 * e * cos(E) - 16)
              * cos(2 * g + 2 * h)
          + (4 / sqrt(1 - pow(e, 2)))
                * (15 * (2 * pow(e, 2) - 1) * cos(E) + 3 * e * (1 - 3 * pow(e, 2)) * cos(2 * E))
                * sin(2 * g + 2 * h);

    // B_M_E
    double B_M_E_0 = 3 * e * (3 * (pow(e, 2) + 4) * cos(E) + e * (e * cos(3 * E) - 6 * cos(2 * E)));
    double B_M_E_1 = -6 * e * (3 * (pow(e, 2) + 4) * cos(E) + e * (e * cos(3 * E) - 6 * cos(2 * E)))
                         * cos(2 * h)
                     + 6
                           * (e * (pow(e, 2) - 2) * cos(3 * E)
                              + (pow(e, 2) + 2) * (2 * cos(2 * E) - 5 * e * cos(E)))
                           * cos(2 * g)
                     + 12 * sqrt(1 - pow(e, 2))
                           * (-2 * (pow(e, 2) + 1) * sin(2 * E) + 5 * e * sin(E) + e * sin(3 * E))
                           * sin(2 * g);
    double B_M_E_2 = 6
                         * ((pow(e, 2) + 2) * (2 * cos(2 * E) - 5 * e * cos(E))
                            + e * (pow(e, 2) - 2) * cos(3 * E))
                         * cos(2 * g - 2 * h)
                     + 12 * sqrt(1 - pow(e, 2))
                           * (-2 * (pow(e, 2) + 1) * sin(2 * E) + 5 * e * sin(E) + e * sin(3 * E))
                           * sin(2 * g - 2 * h);
    double B_M_E_3 = 6
                         * ((pow(e, 2) + 2) * (2 * cos(2 * E) - 5 * e * cos(E))
                            + e * (pow(e, 2) - 2) * cos(3 * E))
                         * cos(2 * g + 2 * h)
                     + 12 * sqrt(1 - pow(e, 2))
                           * (-2 * (pow(e, 2) + 1) * sin(2 * E) + 5 * e * sin(E) + e * sin(3 * E))
                           * sin(2 * g + 2 * h);
    //  double B_M_E_4
    //      = 4
    //        * (9 * e * (4 + pow(e, 2)) * sin(E) - 9 * pow(e, 2) * sin(2 * E) + pow(e, 3) * sin(3 *
    //        E))
    //        * sin(2 * h);

    // B_M_g
    double B_M_g_0 = 0;
    double B_M_g_1 = 4
                         * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                            + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                            - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                         * cos(2 * g)
                     - 4
                           * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                              - e * (2 - pow(e, 2)) * sin(3 * E))
                           * sin(2 * g);
    double B_M_g_2 = -4
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E))
                         * sin(2 * g - 2 * h)
                     + 4
                           * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                              + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                              - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                           * cos(2 * g - 2 * h);
    double B_M_g_3 = -4
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E))
                         * sin(2 * g + 2 * h)
                     + 4
                           * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                              + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                              - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                           * cos(2 * g + 2 * h);

    double B_M_h_0 = 0;
    double B_M_h_1
        = 4
          * (9 * e * (4 + pow(e, 2)) * sin(E) - 9 * pow(e, 2) * sin(2 * E) + pow(e, 3) * sin(3 * E))
          * sin(2 * h);
    double B_M_h_2 = 4
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E))
                         * sin(2 * g - 2 * h)
                     - 4
                           * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                              + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                              - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                           * cos(2 * g - 2 * h);
    double B_M_h_3 = -4
                         * (-15 * e * (2 + pow(e, 2)) * sin(E) + 3 * (2 + pow(e, 2)) * sin(2 * E)
                            - e * (2 - pow(e, 2)) * sin(3 * E))
                         * sin(2 * g + 2 * h)
                     + 4
                           * (-30 * e * sqrt(1 - pow(e, 2)) * cos(E)
                              + 6 * sqrt(1 - pow(e, 2)) * (1 + pow(e, 2)) * cos(2 * E)
                              - 2 * e * sqrt(1 - pow(e, 2)) * cos(3 * E))
                           * cos(2 * g + 2 * h);

    double L_sp2 = (J2_MOON * pow(R_MOON, 2) * sqrt(GM_MOON))
                       / (4 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 3 / 2.0))
                       * (1 - 3 * pow(cos(i), 2)
                          + (2 * pow(1 + e * cos(f), 3) / pow(1 - pow(e, 2), 3 / 2.0))
                                * (1 - 3 * pow(sin(i), 2) * pow(sin(f + g), 2)))
                   + (3 * C22_MOON * pow(R_MOON, 2) * sqrt(GM_MOON))
                         / (2 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 3 / 2.0))
                         * (2 * pow(1 + e * cos(f), 3) / pow(1 - pow(e, 2), 3 / 2.0))
                         * (2 * pow(sin(h) * cos(i) * sin(f + g) - cos(h) * cos(f + g), 2)
                            + pow(sin(i), 2) * pow(sin(f + g), 2) - cos(2 * h) * pow(sin(i), 2))
                   - (k * pow(n3, 2) * pow(a, 7 / 2.0)) / (16 * sqrt(GM_MOON))
                         * (8 * pow(1 - pow(e, 2), 2) / pow(1 + e * cos(f), 2))
                         * (3 * pow(sin(h) * cos(i) * sin(f + g) - cos(h) * cos(f + g), 2) - 1
                            - 15 * pow(e, 2)
                                  * (2 * pow(cos(i / 2), 4) * cos(2 * g + 2 * h)
                                     + 2 * pow(sin(i / 2), 4) * cos(2 * g - 2 * h)
                                     + cos(2 * g) * pow(sin(i), 2))
                            - (3 * pow(e, 2) + 2)
                                  * (3 * cos(2 * h) * pow(sin(i), 2) + 3 * pow(cos(i), 2) - 1));

    double Gsp2
        = ((J2_MOON * sqrt(GM_MOON) * pow(R_MOON, 2))
           / (4 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 3 / 2.0)))
              * (pow(sin(i), 2)
                 * (3 * e * cos(f + 2 * g) + 3 * cos(2 * f + 2 * g) + e * cos(3 * f + 2 * g)))
          +

          ((C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2))
           / (4 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 3 / 2.0)))
              * ((1 - cos(i)) * (1 - cos(i))
                     * (3 * e * cos(f + 2 * g - 2 * h) + 3 * cos(2 * f + 2 * g - 2 * h)
                        + e * cos(3 * f + 2 * g - 2 * h))
                 + (1 + cos(i)) * (1 + cos(i))
                       * (3 * e * cos(f + 2 * g + 2 * h) + 3 * cos(2 * f + 2 * g + 2 * h)
                          + e * cos(3 * f + 2 * g + 2 * h)))
          +

          ((-15 * k * pow(n3, 2) * pow(a, 7 / 2.0) * pow(e, 2)) / (16 * sqrt(GM_MOON)))
              * ((2 * pow(sin(i), 2) * sin(2 * g))
                 + ((1 - cos(i)) * (1 - cos(i)) * sin(2 * g - 2 * h))
                 + ((1 + cos(i)) * (1 + cos(i)) * sin(2 * g + 2 * h)))
              * (E - l)
          +

          ((k * pow(n3, 2) * pow(a, 7 / 2.0)) / (384 * sqrt(GM_MOON)))
              * (4 * (1 - 3 * pow(cos(i), 2)) * B_M_g_0 + 6 * pow(sin(i), 2) * B_M_g_1
                 + 3 * pow(1 - cos(i), 2) * B_M_g_2 + 3 * pow(1 + cos(i), 2) * B_M_g_3);

    double Hsp2 = ((C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2))
                   / (4 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 3 / 2.0)))
                      * (6 * pow(sin(i), 2)
                             * (-2 * (f - l) * sin(2 * h) - e * cos(f - 2 * h) + e * cos(f + 2 * h))
                         - pow(1 - cos(i), 2)
                               * (3 * e * cos(f + 2 * g - 2 * h) + 3 * cos(2 * f + 2 * g - 2 * h)
                                  + e * cos(3 * f + 2 * g - 2 * h))
                         + pow(1 + cos(i), 2)
                               * (3 * e * cos(f + 2 * g + 2 * h) + 3 * cos(2 * f + 2 * g + 2 * h)
                                  + e * cos(3 * f + 2 * g + 2 * h)))
                  + ((3 * k * pow(n3, 2) * pow(a, 7 / 2.0)) / (16 * sqrt(GM_MOON)))
                        * (-2 * pow(sin(i), 2) * (2 + 3 * pow(e, 2)) * sin(2 * h)
                           + 5 * pow(e, 2) * pow(1 - cos(i), 2) * sin(2 * g - 2 * h)
                           - 5 * pow(e, 2) * pow(1 + cos(i), 2) * sin(2 * g + 2 * h))
                        * (E - l)
                  + ((k * pow(n3, 2) * pow(a, 7 / 2.0)) / (384 * sqrt(GM_MOON)))
                        * (4 * (1 - 3 * pow(cos(i), 2)) * B_M_h_0 + 6 * pow(sin(i), 2) * B_M_h_1
                           + 3 * pow(1 - cos(i), 2) * B_M_h_2 + 3 * pow(1 + cos(i), 2) * B_M_h_3);

    double lsp2
        = (-(J2_MOON * pow(R_MOON, 2)) / (32 * pow(a, 2) * e * pow(1 - pow(e, 2), 3 / 2.0)))
              * (-4 * (1 - 3 * pow(cos(i), 2))
                     * (pow(e, 2) * cos(2 * f) - pow(e, 2) + 6 * e * cos(f) + 6) * sin(f)
                 + pow(sin(i), 2)
                       * (-18 * e * sin(2 * g) - 3 * (4 + 5 * pow(e, 2)) * sin(f + 2 * g)
                          + 3 * pow(e, 2) * sin(f - 2 * g) + (28 - pow(e, 2)) * sin(3 * f + 2 * g)
                          + 18 * e * sin(4 * f + 2 * g) + 3 * pow(e, 2) * sin(5 * f + 2 * g)))
          + (-(C22_MOON * pow(R_MOON, 2)) / (32 * pow(a, 2) * e * pow(1 - pow(e, 2), 3 / 2.0)))
                * (24 * pow(sin(i), 2) * sin(f)
                       * (pow(e, 2) * cos(2 * f) - pow(e, 2) + 6 * e * cos(f) + 6) * cos(2 * h)
                   + pow(1 - cos(i), 2)
                         * (-3 * (5 * pow(e, 2) + 4) * sin(f + 2 * g - 2 * h)
                            - (pow(e, 2) - 28) * sin(3 * f + 2 * g - 2 * h)
                            + 6 * e * sin(f) * (2 * e * cos(2 * f) + e + 12 * cos(f))
                                  * cos(2 * f + 2 * g - 2 * h))
                   + pow(1 + cos(i), 2)
                         * (pow(e, 2)
                                * (-15 * sin(f + 2 * g + 2 * h) - sin(3 * f + 2 * g + 2 * h)
                                   + 3 * sin(5 * f + 2 * g + 2 * h) + 3 * sin(f - 2 * g - 2 * h))
                            + 18 * e * sin(4 * f + 2 * g + 2 * h) - 18 * e * sin(2 * g + 2 * h)
                            - 12 * sin(f + 2 * g + 2 * h) + 28 * sin(3 * f + 2 * g + 2 * h)))
          + (-(k * pow(n3, 2) * pow(a, 3)) / (384 * (1 - e * sin(E))))
                * (84
                       * (-2 * (1 - 3 * pow(cos(i), 2)) * (2 + 3 * pow(e, 2))
                          + 6 * pow(sin(i), 2) * (2 + 3 * pow(e, 2)) * cos(2 * h)
                          + 30 * pow(e, 2) * pow(sin(i), 2) * cos(2 * g)
                          + 15 * pow(e, 2) * pow(1 - cos(i), 2) * cos(2 * g - 2 * h)
                          + 15 * pow(e, 2) * pow(1 + cos(i), 2) * cos(2 * g + 2 * h))
                       * (E - l)
                   + 28 * (1 - 3 * pow(cos(i), 2)) * B_M_1 + 42 * pow(sin(i), 2) * B_M_2
                   + 21 * pow(1 - cos(i), 2) * B_M_3 + 21 * pow(1 + cos(i), 2) * B_M_4
                   - 72 * (pow(e, 2) - 1)
                         * (-2 * (1 - 3 * pow(cos(i), 2)) + 6 * pow(sin(i), 2) * cos(2 * h)
                            + 10 * pow(sin(i), 2) * cos(2 * g)
                            + 5 * pow(1 - cos(i), 2) * cos(2 * g - 2 * h)
                            + 5 * pow(1 + cos(i), 2) * cos(2 * g + 2 * h))
                   + (1 - pow(e, 2)) / e
                         * (4 * (1 - 3 * pow(cos(i), 2)) * B_M_e_0 + 6 * pow(sin(i), 2) * B_M_e_1
                            + 3 * pow(1 - cos(i), 2) * B_M_e_2 + 3 * pow(1 + cos(i), 2) * B_M_e_3)
                   - 24 * (1 - 3 * pow(cos(i), 2)) * (2 + 3 * pow(e, 2))
                   + 72 * pow(sin(i), 2) * (2 + 3 * pow(e, 2)) * cos(2 * h)
                   + 360 * pow(e, 2) * pow(sin(i), 2) * cos(2 * g)
                   + 180 * pow(e, 2) * pow(1 - cos(i), 2) * cos(2 * g - 2 * h)
                   + 180 * pow(e, 2) * pow(1 + cos(i), 2) * cos(2 * g + 2 * h)
                   + 4 * (1 - 3 * pow(cos(i), 2)) * B_M_E_0 + 6 * pow(sin(i), 2) * B_M_E_1
                   + 3 * pow(1 - cos(i), 2) * B_M_E_2 + 3 * pow(1 + cos(i), 2) * B_M_E_3)
                * (sin(E) / (1 - e * sin(E)));

    double gsp2
        = ((J2_MOON * pow(R_MOON, 2)) / (32 * pow(a, 2) * e * pow(1 - pow(e, 2), 2)))
              * ((pow(sin(i), 2)
                  * (-18 * e * sin(2 * g) - 3 * (4 - 7 * pow(e, 2)) * sin(f + 2 * g)
                     + 3 * pow(e, 2) * sin(f - 2 * g) + 36 * e * sin(2 * f + 2 * g)
                     + (28 + 11 * pow(e, 2)) * sin(3 * f + 2 * g) + 18 * e * sin(4 * f + 2 * g)
                     + 3 * pow(e, 2) * sin(5 * f + 2 * g)))
                 + ((3 * cos(2 * i) + 1)
                    * (3 * (3 * pow(e, 2) + 4) * sin(f)
                       + e * (e * sin(3 * f) + 6 * (sin(2 * f) + 2 * f - 2 * l))))
                 + (-8 * e
                    * (e * (3 * sin(f + 2 * g) + sin(3 * f + 2 * g) - 6 * sin(f)) - 6 * f + 6 * l
                       + 3 * sin(2 * f + 2 * g))
                    * pow(cos(i), 2)))
          + ((C22_MOON * pow(R_MOON, 2)) / (32 * pow(a, 2) * e * pow(1 - pow(e, 2), 2)))
                * ((24 * pow(sin(i), 2)
                    * (sin(f) * (pow(e, 2) * cos(2 * f) + 5 * pow(e, 2) + 6 * e * cos(f) + 6)
                       + 6 * e * (f - l))
                    * cos(2 * h))
                   + (-8 * e * pow(cos(i), 2)
                      * (6 * cos(i) * B_M_1 + (1 - cos(i)) * B_M_2 - (1 + cos(i)) * B_M_3))
                   + (pow(1 - cos(i), 2)
                      * (-18 * e * sin(2 * g - 2 * h)
                         - 3 * (4 - 7 * pow(e, 2)) * sin(f + 2 * g - 2 * h)
                         + 3 * pow(e, 2) * sin(f - 2 * g + 2 * h)
                         + 36 * e * sin(2 * f + 2 * g - 2 * h)
                         + (28 + 11 * pow(e, 2)) * sin(3 * f + 2 * g - 2 * h)
                         + 18 * e * sin(4 * f + 2 * g - 2 * h)
                         + 3 * pow(e, 2) * sin(5 * f + 2 * g - 2 * h)))
                   + (pow(1 + cos(i), 2)
                      * (-18 * e * sin(2 * g + 2 * h)
                         - 3 * (4 - 7 * pow(e, 2)) * sin(f + 2 * g + 2 * h)
                         + 3 * pow(e, 2) * sin(f - 2 * g - 2 * h)
                         + 36 * e * sin(2 * f + 2 * g + 2 * h)
                         + (28 + 11 * pow(e, 2)) * sin(3 * f + 2 * g + 2 * h)
                         + 18 * e * sin(4 * f + 2 * g + 2 * h)
                         + 3 * pow(e, 2) * sin(5 * f + 2 * g + 2 * h))))
          + ((k * pow(n3, 2) * pow(a, 3)) / (64 * GM_MOON * sqrt(1 - pow(e, 2))))
                * (((E - l)
                        * ((24 * pow(cos(i), 2)
                            * (10 * pow(e, 2) * cos(2 * g) * pow(sin(h), 2)
                               + (3 * pow(e, 2) + 2) * (cos(h) - 1)))
                           + (120 * pow(e, 2) * cos(i) * sin(2 * g) * sin(2 * h)
                              + 12 * (pow(e, 2) - 1)
                                    * (2 * (3 - 5 * cos(2 * g)) * pow(sin(h), 2) * cos(2 * i)
                                       - 20 * sin(2 * g) * sin(2 * h) * cos(i)
                                       + (5 * cos(2 * g) + 1) * (3 * cos(2 * h) + 1))))
                    + 4 * pow(cos(i), 2) * B_M_1 + 2 * pow(cos(i), 2) * B_M_2
                    + (1 - cos(i)) * cos(i) * B_M_3 - (1 + cos(i)) * cos(i) * B_M_4)
                   + (((pow(e, 2) - 1) / (6 * e))
                      * (4 * (1 - 3 * pow(cos(i), 2)) * (B_M_E_0) + 6 * pow(sin(i), 2) * (B_M_E_1)
                         + 3 * pow(1 - cos(i), 2) * (B_M_E_2) + 3 * pow(1 + cos(i), 2) * (B_M_E_3)
                         + (-24 * (1 - 3 * pow(cos(i), 2)) * (2 + 3 * pow(e, 2))
                            + 72 * pow(sin(i), 2) * (2 + 3 * pow(e, 2)) * cos(2 * h)
                            + 360 * pow(e, 2) * pow(sin(i), 2) * cos(2 * g)
                            + 180 * pow(e, 2) * pow(1 - cos(i), 2) * cos(2 * g - 2 * h)
                            + 180 * pow(e, 2) * pow(1 + cos(i), 2) * cos(2 * g + 2 * h)
                            + 4 * (1 - 3 * pow(cos(i), 2)) * (B_M_E_0)
                            + 6 * pow(sin(i), 2) * (B_M_E_1) + 3 * pow(1 - cos(i), 2) * (B_M_E_2)
                            + 3 * pow(1 + cos(i), 2) * (B_M_E_3))
                               * (sin(E) / (1 - e * sin(E))))));

    double hsp2 = ((J2_MOON * pow(R_MOON, 2) * cos(i)) / (4 * pow(a, 2) * pow(1 - pow(e, 2), 2))
                   * (-6 * B_2_1 + B_2_2))
                  + ((C22_MOON * pow(R_MOON, 2)) / (4 * pow(a, 2) * pow(1 - pow(e, 2), 2))
                     * (6 * cos(i) * B_22_1 + (1 - cos(i)) * B_22_2 - (1 + cos(i)) * B_22_3))
                  + ((k * pow(n3, 2) * pow(a, 3)) / (64 * GM_MOON * sqrt(1 - pow(e, 2))))
                        * ((-2 * cos(i) * (2 + 3 * pow(e, 2))
                            + 2 * cos(i) * (2 + 3 * pow(e, 2)) * cos(2 * h)
                            + 10 * pow(e, 2) * cos(i) * cos(2 * g)
                            + 5 * pow(e, 2) * (1 - cos(i)) * cos(2 * g - 2 * h)
                            - 5 * pow(e, 2) * (1 + cos(i)) * cos(2 * g + 2 * h))
                               * (E - l)
                           + (4 * cos(i) * B_M_1 + 2 * cos(i) * B_M_2 + (1 - cos(i)) * B_M_3
                              - (1 + cos(i)) * B_M_4));

    std::array<double, 6> ret = {L_sp2, Gsp2, Hsp2, lsp2, gsp2, hsp2};
    return ret;
  }

  std::array<double, 6> ComputeFirstOrderMediumPeriod(Vec6 &coe, Vec6 &doe) {
    double n3 = 2.66e-6;
    //  double nM = 2.66e-6;
    //  double J2 = 2.03e-4;
    double k = 0.98785;

    double a = (double)coe(0);
    double e = (double)coe(1);
    double i = (double)coe(2);
    //  double O = (double)coe(3);
    //  double w = (double)coe(4);
    //  double M = (double)coe(5);

    //  double l = (double)doe(0);
    double g = (double)doe(1);
    double h = (double)doe(2);
    //  double L = (double)doe(3);
    //  double G = (double)doe(4);
    //  double H = (double)doe(5);

    // The first-order terms of the medium-period variation are
    double L_mp1 = 0;

    double G_mp1
        = -15 * k * n3 / 32 * pow(a, 2) * pow(e, 2)
          * (4 * cos(i) * cos(2 * g) * cos(2 * h) - (cos(2 * i) + 3) * sin(2 * g) * sin(2 * h));

    double H_mp1 = -3 * k * n3 * pow(a, 2) / 32
                       * ((5 * pow(e, 2) * (cos(2 * i) + 3) * cos(2 * g)
                           + 2 * (3 * pow(e, 2) + 2) * pow(sin(i), 2))
                              * cos(2 * h)
                          - 20 * pow(e, 2) * cos(i) * sin(2 * g) * sin(2 * h))
                   - 3 * C22_MOON * GM_MOON * pow(R_MOON, 2)
                         / (2 * n3 * pow(a, 3) * pow(1 - pow(e, 2), 3 / 2.0)) * pow(sin(i), 2)
                         * cos(2 * h);

    double l_mp1 = -9 * C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * pow(sin(i), 2) * sin(2 * h)
                       / (4 * n3 * pow(a, 7 / 2.0) * pow(1 - pow(e, 2), 3 / 2.0))
                   + 3 * k * n3 * pow(a, 3 / 2.0) / (32 * sqrt(GM_MOON))
                         * (sin(2 * h)
                                * (5 * (pow(e, 2) + 1) * cos(2 * g) * (cos(2 * i) + 3)
                                   + 2 * (3 * pow(e, 2) + 7) * pow(sin(i), 2))
                            + 20 * (pow(e, 2) + 1) * sin(2 * g) * cos(2 * h) * cos(i));

    double g_mp1 = 3 * C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * (5 * cos(2 * i) - 1) * sin(2 * h)
                       / (8 * n3 * pow(a, 7 / 2.0) * pow(1 - pow(e, 2), 2))
                   + 3 * k * n3 * pow(a, 3 / 2.0) / (64 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                         * (20 * (pow(e, 2) - 2) * sin(2 * g) * cos(2 * h) * cos(i)
                            + (10 * (2 * pow(e, 2) - 3) * cos(2 * g) + 12 * pow(e, 2)
                               + 20 * pow(sin(g), 2) * cos(2 * i) - 2)
                                  * sin(2 * h));

    double h_mp1
        = -3 * C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * cos(i) * sin(2 * h)
              / (2 * n3 * pow(a, 7 / 2.0) * pow(1 - pow(e, 2), 2))
          + 3 * k * n3 * pow(a, 3 / 2.0) / (16 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                * (5 * pow(e, 2) * sin(2 * g) * cos(2 * h)
                   + sin(2 * h) * cos(i) * (5 * pow(e, 2) * cos(2 * g) - 3 * pow(e, 2) - 2));

    std::array<double, 6> ret = {l_mp1, g_mp1, h_mp1, L_mp1, G_mp1, H_mp1};
    return ret;
  }

  std::array<double, 6> ComputeSecondOrderMediumPeriod(Vec6 &coe, Vec6 &doe) {
    double n3 = 2.66e-6;
    //  double nM = 2.66e-6;
    //  double J2 = 2.03e-4;
    double k = 0.98785;

    double a = (double)coe(0);
    double e = (double)coe(1);
    double i = (double)coe(2);
    //  double O = (double)coe(3);
    //  double w = (double)coe(4);
    //  double M = (double)coe(5);

    //  double l = (double)doe(0);
    double g = (double)doe(1);
    double h = (double)doe(2);
    //  double L = (double)doe(3);
    //  double G = (double)doe(4);
    //  double H = (double)doe(5);

    // The second-order terms of the medium-period variation are
    double L_mp2 = 0;

    double G_mp2
        = -45 * k * J2_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * pow(e, 2)
              / (512 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 2))
              * (cos(2 * g) * cos(2 * h) * (20 * cos(2 * i) + 5 * cos(4 * i) + 7)
                 - 16 * sin(2 * g) * sin(2 * h) * (cos(i) + cos(3 * i)))
          + 45 * k * C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * pow(e, 2) * sin(i)
                / (256 * pow(a, 3 / 2.0) * pow(1 - pow(e, 2), 2))
                * (2 * sin(2 * g) * (4 * sin(2 * h) - 5 * sin(4 * h)) * sin(2 * i)
                   + cos(2 * g) * sin(i)
                         * (5 * (-4 * cos(2 * h) + cos(4 * h) - 2) * cos(2 * i) + 4 * cos(2 * h)
                            + 15 * cos(4 * h) - 22))
          + 45 * pow(k, 2) * pow(n3, 2) * pow(a, 7 / 2.0) * pow(e, 2)
                / (4096 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                * (-5 * (-40 * pow(e, 2) + 3 * cos(4 * g) - 3 * cos(4 * h) + 28 * cos(2 * i) + 35)
                   + 4 * cos(2 * i)
                         * (5 * (6 * pow(e, 2) - 2 * pow(cos(g), 2) * cos(4 * h) + cos(4 * g))
                            - 32 * (pow(e, 2) - 1) * sin(2 * g) * sin(2 * h))
                   + cos(2 * g)
                         * (8 * (pow(e, 2) - 1) * (16 * cos(2 * h) - 3) * cos(2 * i)
                            + 24 * pow(e, 2) + 20 * pow(sin(g), 2) * cos(4 * i) + 25 * cos(4 * h)
                            - 34)
                   + 10
                         * (pow(sin(g), 2) * cos(4 * h) * cos(4 * i)
                            - 8 * sin(2 * g) * sin(4 * h) * pow(sin(i), 2) * cos(i)))
          - pow(e, 2) / (1024 * sqrt(GM_MOON) * pow(a, 3 / 2.0))
                * (C22_MOON * GM_MOON * pow(R_MOON, 2) * (105 * pow(e, 4) + 80 * pow(e, 2) + 48)
                   - 480 * k * pow(n3, 2) * pow(a, 5) * sqrt(1 - pow(e, 2)))
                * ((cos(2 * g) * cos(2 * h) * (cos(2 * i) + 3))
                   - (4 * sin(2 * g) * sin(2 * h) * cos(i)));

    double H_mp2
        = (9 * J2_MOON * C22_MOON * pow(GM_MOON, 1.5) * pow(R_MOON, 4) * pow(sin(i), 2) * cos(i))
              / (4 * pow(n3, 2) * pow(a, 6.5) * pow(1 - pow(e, 2), 3.5)) * cos(2 * h)
          - (9 * J2_MOON * k * pow(R_MOON, 2) * sqrt(GM_MOON))
                / (512 * sqrt(a) * pow(1 - pow(e, 2), 2))
                * (80 * pow(e, 2) * cos(2 * g) * cos(2 * h) * (cos(i) + cos(3 * 2 * sin(i)))
                   - 5 * pow(e, 2) * sin(2 * g) * sin(2 * h)
                         * (20 * cos(2 * i) + 5 * cos(4 * 2 * cos(i)) + 7)
                   - 16 * (3 * pow(e, 2) + 2) * cos(2 * h) * pow(sin(i), 2) * cos(i))
          - (9 * pow(C22_MOON, 2) * pow(GM_MOON, 1.5) * pow(R_MOON, 4) * pow(sin(i), 2) * cos(i))
                / (4 * pow(n3, 2) * pow(a, 6.5) * pow(1 - pow(e, 2), 3.5))
          - (9 * k * C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * pow(sin(i), 2))
                / (128 * sqrt(a) * pow(1 - pow(e, 2), 2))
                * (20 * pow(e, 2) * cos(2 * g) * (2 * cos(2 * h) + 3) * cos(i)
                   - 50 * pow(e, 2) * sin(2 * g) * sin(2 * h) * pow(cos(i), 2)
                   + 5 * pow(e, 2) * sin(2 * g) * sin(2 * h) * (7 - 5 * cos(2 * i))
                   + 16 * (3 * pow(e, 2) + 2) * pow(sin(2 * h), 2) * cos(i))
          + (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3.5) * sqrt(1 - pow(e, 2)))
                / (1024 * sqrt(GM_MOON))
                * (-160 * pow(e, 2) * cos(2 * i) * cos(2 * (cos(2 * g) + cos(2 * h)))
                   - 2 * cos(i)
                         * (15 * pow(e, 2) * cos(2 * g) + (34 * pow(e, 2) - 4) * cos(2 * h)
                            + 183 * pow(e, 2) + 2)
                   + 2 * cos(3 * 2 * sin(i))
                         * (15 * pow(e, 2) * cos(2 * g)
                            + (17 * pow(e, 2) - 2) * (2 * cos(2 * h) - 1)))
          - (pow(e, 2)) / (1024 * sqrt(GM_MOON) * sqrt(a))
                * (C22_MOON * GM_MOON * pow(R_MOON, 2) * (105 * pow(e, 4) + 80 * pow(e, 2) + 48)
                   - 480 * k * pow(n3, 2) * pow(a, 5) * sqrt(1 - pow(e, 2)))
                * (4 * cos(2 * g) * cos(2 * h) * cos(i)
                   - sin(2 * g) * sin(2 * h) * (cos(2 * i) + 3));

    // TODO check parenthesis
    double l_mp2
        = (27 * J2_MOON * C22_MOON * GM_MOON * pow(R_MOON, 4) * pow(sin(i), 2) * cos(i))
              / (4 * pow(n3, 2) * pow(a, 7) * pow(1 - pow(e, 2), 7 / 2)) * sin(2 * h)
          - (9 * J2_MOON * k * pow(R_MOON, 2)) / (1024 * pow(a, 2) * pow(1 - pow(e, 2), 2))
                * (5 * (pow(e, 2) - 2) * sin(2 * g) * cos(2 * h)
                       * (20 * cos(2 * i) + 5 * cos(4 * i) + 7)
                   + 8 * sin(2 * h)
                         * (10 * (pow(e, 2) - 2) * cos(2 * g) * (cos(i) + cos(3 * i))
                            + (8 - 3 * pow(e, 2)) * sin(i) * sin(2 * i)))
          - (9 * k * C22_MOON * pow(R_MOON, 2) * pow(sin(i), 2))
                / (512 * pow(a, 2) * pow(1 - pow(e, 2), 2))
                * (5 * sin(2 * g)
                   * (4 * (pow(e, 2) - 2) * cos(2 * h) * (5 * cos(2 * i) - 1)
                      + 5 * cos(4 * h) * ((pow(e, 2) - 2) * cos(2 * i) + 15 * pow(e, 2) + 6)
                      - 4
                            * (10 * (pow(e, 2) + 1) * pow(sin(2 * h), 2) * cos(2 * i)
                               + 5 * pow(e, 2) + 11)
                      + 4 * sin(2 * h) * cos(i)
                            * (125 * pow(e, 2) * cos(2 * g + 2 * h)
                               + 25 * (5 * pow(e, 2) + 2) * cos(2 * g - 2 * h)
                               + 20 * (pow(e, 2) - 2) * cos(2 * g) - 12 * pow(e, 2)
                               + 50 * cos(2 * g + 2 * h) + 32)))
          +  // Added parenthesis here
          (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3)) / (8192 * GM_MOON * pow(1 - pow(e, 2), 3 / 2))
              * ((pow(e, 2) - 1) * (pow(e, 2) - 1)
                     * (8
                            * (160 * (2 * pow(e, 2) + 1) * cos(2 * i) * sin(2 * g + 2 * h)
                               + 8 * (34 * pow(e, 2) + 11) * sin(2 * h) * sin(i) * sin(2 * i)
                               - 25 * sin(4 * h)
                                     * (cos(4 * g) * (7 * cos(i) + cos(3 * i))
                                        + 4 * cos(2 * g) * pow(sin(i), 2) * cos(i)))
                        + 25 * cos(4 * h)
                              * (2 * sin(2 * g) * (4 * cos(2 * i) - 5)
                                 - 7 * sin(4 * g) * (4 * cos(2 * i) + 5)))
                 + 40 * sin(g) * cos(g)
                       * (5 * cos(2 * g) * (4 * (pow(e, 4) - 1) * cos(2 * i) - 3 * pow(e, 4))
                          - 33 * pow(e, 4)
                          + 5 * (pow(e, 2) - 1) * pow(sin(g), 2) * cos(4 * i)
                                * ((pow(e, 2) - 1) * cos(4 * h) + 2 * (pow(e, 2) + 1))
                          + 16 * pow(e, 2) + 4 * (7 * pow(e, 4) - 4 * pow(e, 2) - 3) * cos(2 * i))
                 + 10 * (34 * sin(2 * g) + 15 * sin(4 * g)))
          - (1 / (2048 * pow(a, 2) * GM_MOON))
                * (960 * k * pow(n3, 2) * pow(a, 5) * (2 * pow(e, 2) + 1) * sqrt(1 - pow(e, 2))
                   + C22_MOON * GM_MOON * pow(R_MOON, 2)
                         * (945 * pow(e, 6) - 70 * pow(e, 4) - 80 * pow(e, 2) - 96))
                * (sin(2 * g) * cos(2 * h) * (cos(2 * i) + 3)
                   + 4 * cos(2 * g) * sin(2 * h) * cos(i));

    double g_mp2
        = -(9 * J2_MOON * C22_MOON * GM_MOON * pow(R_MOON, 4) * cos(i) * (5 * cos(2 * i) - 3))
              / (8 * pow(n3, 2) * pow(a, 7) * pow(1 - pow(e, 2), 4)) * sin(2 * h)
          - (9 * k * J2_MOON * pow(R_MOON, 2)) / (512 * pow(a, 2) * pow(1 - pow(e, 2), 5 / 2))
                * (5 * sin(2 * g) * cos(2 * h)
                       * (20 * (3 * pow(e, 2) + 1) * cos(2 * i)
                          + 5 * (3 * pow(e, 2) + 1) * cos(4 * i) + 37 * pow(e, 2) + 7)
                   + 4 * sin(2 * h) * cos(i)
                         * (5 * cos(2 * i)  // Break
                                * (4 * (5 * pow(e, 2) + 2) * cos(2 * g) + 3 * pow(e, 2) + 4)
                            + 40 * pow(e, 2) * cos(2 * g) - 3 * pow(e, 2) - 12))
          - (9 * pow(C22_MOON, 2) * GM_MOON * pow(R_MOON, 4) * sin(4 * h) * cos(i))
                / (8 * pow(n3, 2) * pow(a, 7) * pow(1 - pow(e, 2), 4))
          + (9 * k * C22_MOON * pow(R_MOON, 2)) / (1024 * pow(a, 2) * pow(1 - pow(e, 2), 5 / 2))
                * (-8
                       * (10 * pow(e, 2) * sin(2 * g) * pow(cos(2 * h), 2) * (1 - 5 * cos(2 * i))
                          + 10 * sin(2 * g) * pow(sin(2 * h), 2) * pow(sin(i), 2)
                                * (-6 * pow(e, 2) + 5 * cos(2 * i) + 11)
                          + sin(4 * h) * cos(i)
                                * (10 * cos(2 * g)
                                       * ((4 - 5 * pow(e, 2)) * cos(2 * i) + 3 * pow(e, 2) - 4)
                                   + 5 * (3 * pow(e, 2) + 4) * cos(2 * i) + 9 * pow(e, 2) - 4))
                   - 4 * cos(i)
                         * (16 * pow(sin(h), 3) * cos(h)
                                * (5 * (3 * pow(e, 2) + 4) * cos(2 * i) - 3 * (pow(e, 2) + 4))
                            - 5 * cos(2 * g) * (4 * sin(2 * h) + 3 * sin(4 * h))
                                  * ((5 * pow(e, 2) + 2) * cos(2 * i) - pow(e, 2) - 2))
                   + 5 * sin(2 * g) * cos(4 * h)
                         * (4 * (7 * pow(e, 2) + 1) * cos(2 * i)
                            + 5 * (3 * pow(e, 2) + 1) * cos(4 * i) + 5 * pow(e, 2) - 9)
                   + 20 * sin(2 * g) * cos(2 * h)
                         * (-4 * (pow(e, 2) + 3) * cos(2 * i) + 5 * pow(e, 2)
                            + 5 * (3 * pow(e, 2) + 1) * cos(4 * i) + 7))
          + (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3)) / (16384 * GM_MOON * (1 - pow(e, 2)))
                * (10 * sin(2 * g) * cos(4 * h)
                       * (-120 * pow(e, 4) + 5 * (3 * pow(e, 2) + 2) * cos(4 * i) + 133 * pow(e, 2)
                          + 4 * (6 * pow(e, 4) - 13 * pow(e, 2) + 2) * cos(2 * i) - 18)
                   - 25 * sin(4 * g) * cos(4 * h)
                         * (72 * pow(e, 4) + (3 * pow(e, 2) + 2) * cos(4 * i) - 135 * pow(e, 2)
                            + 4 * (6 * pow(e, 4) - 15 * pow(e, 2) + 14) * cos(2 * i) + 70)
                   + 64 * (pow(e, 2) - 1) * sin(2 * h)
                         * (20 * cos(2 * g) * (2 * pow(e, 2) - (pow(e, 2) - 2) * cos(2 * i))
                            - 2 * cos(i) * (34 * pow(e, 2) + 15 * cos(2 * i) - 19))
                   + 320 * (pow(e, 2) - 1) * cos(2 * h)
                         * (sin(2 * g) * (8 * pow(e, 2) - 4 * (pow(e, 2) - 2) * cos(2 * i))
                            - 3 * cos(2 * g) * sin(2 * h) * cos(i)
                                  * (2 * pow(e, 2) + cos(2 * i) - 1))
                   - 16
                         * (2
                                * (2 * sin(2 * h) * cos(i)
                                       * (5 * pow(e, 2) * cos(2 * g) - 3 * pow(e, 2) - 2)
                                   + 10 * pow(e, 2) * sin(2 * g) * cos(2 * h))
                                * (cos(2 * h)
                                       * (5 * cos(2 * g) * (-2 * pow(e, 2) + cos(2 * i) + 3)
                                          - 6 * pow(e, 2) - 5 * cos(2 * i) + 1)
                                   + 10 * (pow(e, 2) - 2) * sin(2 * g) * sin(2 * h) * cos(i))
                            - 5
                                  * (2 * sin(2 * h)
                                         * (5 * (2 * pow(e, 2) - 3) * cos(2 * g) + 6 * pow(e, 2)
                                            + 10 * pow(sin(2 * g), 2) * cos(2 * i) - 1)
                                     + 20 * (pow(e, 2) - 2) * sin(2 * g) * cos(2 * h) * cos(i))
                                  * (2 * (pow(e, 2) - 2) * cos(2 * g) * cos(2 * h) * cos(i)
                                     + sin(2 * g) * sin(2 * h) * (-2 * pow(e, 2) + cos(2 * i) + 3))
                            - 16 * sin(4 * h) * cos(i)
                                  * (50 * cos(4 * g) * (3 * pow(e, 2) - 1) * (pow(e, 2) - 1)
                                     + cos(2 * i)
                                     - 2 * (3 * pow(e, 2) + 2)
                                           * (3 * pow(e, 2) + 5 * cos(2 * i) - 3)))
                   - (15 * k * pow(n3, 2) * pow(a, 3)) / (64 * GM_MOON)
                         * (sin(2 * g) * cos(2 * h)
                                * ((pow(e, 2) - 2) * cos(2 * i) + 7 * pow(e, 2) - 6)
                            + 8 * (pow(e, 2) - 1) * cos(2 * g) * sin(2 * h) * cos(i))
                   - (C22_MOON * pow(R_MOON, 2)) / (1024 * pow(a, 2) * sqrt(1 - pow(e, 2)))
                         * ((1 - pow(e, 2)) * (315 * pow(e, 4) + 160 * pow(e, 2) + 48)
                                * (sin(2 * g) * cos(2 * h) * (cos(2 * i) + 3)
                                   + 4 * cos(2 * g) * sin(2 * h) * cos(i))
                            + (105 * pow(e, 4) + 80 * pow(e, 2) + 48) * pow(e, 2) * (1 / tan(i))
                                  * (2 * cos(2 * g) * sin(2 * h) * sin(i)
                                     + sin(2 * g) * cos(2 * h) * sin(2 * i))));

    double h_mp2
        = (9 * J2_MOON * C22_MOON * GM_MOON * pow(R_MOON, 4) * (3 * cos(2 * i) + 1))
              / (16 * pow(n3, 2) * pow(a, 7) * pow(1 - pow(e, 2), 4)) * sin(2 * h)
          + (9 * k * J2_MOON * pow(R_MOON, 2)) / (128 * pow(a, 2) * pow(1 - pow(e, 2), 5 / 2))
                * (100 * pow(e, 2) * sin(2 * g) * cos(2 * h) * pow(cos(i), 3)
                   + sin(2 * h)
                         * (20 * pow(e, 2) * cos(2 * g) * (3 * cos(2 * i) + 2)
                            + (3 * pow(e, 2) + 2) * (3 * cos(2 * i) + 1)))
          + (9 * pow(C22_MOON, 2) * GM_MOON * pow(R_MOON, 4) * sin(4 * h) * (cos(2 * i) + 3))
                / (32 * pow(n3, 2) * pow(a, 7) * pow(1 - pow(e, 2), 4))
          + (9 * k * C22_MOON * pow(R_MOON, 2)) / (256 * pow(a, 2) * pow(1 - pow(e, 2), 5 / 2))
                * (-2
                       * (5 * pow(e, 2) * sin(2 * g) * cos(i)
                              * (10 * pow(sin(2 * h), 2) * cos(2 * i) + 9 * cos(4 * h) - 1)
                          - 5 * pow(e, 2) * cos(2 * g) * sin(4 * h) * (cos(2 * i) - 5)
                          - 8 * (3 * pow(e, 2) + 2) * sin(4 * h) * pow(cos(i), 2))
                   + (3 * cos(2 * i) + 1)
                         * (16 * (3 * pow(e, 2) + 2) * pow(sin(h), 3) * cos(h)
                            - 5 * pow(e, 2) * cos(2 * g) * (4 * sin(2 * h) + 3 * sin(4 * h)))
                   - 5 * pow(e, 2) * sin(2 * g) * cos(4 * h) * (7 * cos(i) + 5 * cos(3 * i))
                   + 20 * pow(e, 2) * sin(2 * g) * cos(2 * h) * (cos(i) - 5 * cos(3 * i)))
          - (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3)) / (2048 * GM_MOON * (1 - pow(e, 2)))
                * (10 * pow(e, 2) * cos(i)
                       * (2 * sin(2 * g)
                              * (32 * (pow(e, 2) - 1) * cos(2 * h) + 6 * pow(e, 2) + 10 * cos(4 * h)
                                 + 5 * cos(2 * i) - 1)
                          - 5 * sin(4 * g) * (2 * pow(e, 2) * cos(4 * h) + cos(2 * i) - 1)
                          + 64 * (pow(e, 2) - 1) * cos(2 * g) * sin(2 * h))
                   + sin(4 * h)
                         * (-25 * pow(e, 4) * cos(4 * g) * (cos(2 * i) + 3) - 27 * pow(e, 4)
                            + 10 * pow(e, 2) * cos(2 * g)
                                  * ((3 * pow(e, 2) + 7) * cos(2 * i) - 3 * pow(e, 2) + 13)
                            + 14 * pow(e, 2) - (9 * pow(e, 4) + 62 * pow(e, 2) + 4) * cos(2 * i)
                            - 12)
                   - 8 * (17 * pow(e, 4) - 19 * pow(e, 2) + 2) * sin(2 * h) * (3 * cos(2 * i) + 1))
          + (pow(e, 2) / (512 * GM_MOON * pow(a, 2) * sqrt(1 - pow(e, 2))))
                * (C22_MOON * GM_MOON * pow(R_MOON, 2) * (105 * pow(e, 4) + 80 * pow(e, 2) + 48)
                   - 480 * k * pow(n3, 2) * pow(a, 5) * sqrt(1 - pow(e, 2)))
                * (sin(2 * g) * cos(2 * h) * cos(i) + cos(2 * g) * sin(2 * h));

    std::array<double, 6> ret = {l_mp2, g_mp2, h_mp2, L_mp2, G_mp2, H_mp2};
    return ret;
  }

  std::array<double, 6> ComputeCorrectionMediumPeriod(Vec6 &coe, Vec6 &doe) {
    double n3 = 2.66e-6;
    //  double nM = 2.66e-6;
    //  double J2 = 2.03e-4;
    double k = 0.98785;

    double a = (double)coe(0);
    double e = (double)coe(1);
    double i = (double)coe(2);
    //  double O = (double)coe(3);
    //  double w = (double)coe(4);
    //  double M = (double)coe(5);

    //  double l = (double)doe(0);
    double g = (double)doe(1);
    double h = (double)doe(2);
    //  double L = (double)doe(3);
    //  double G = (double)doe(4);

    // Assuming the required variables are already defined
    double L_mp1c = 0;

    double G_mp1c
        = (45 * k * C22_MOON * sqrt(GM_MOON) * pow(R_MOON, 2) * pow(e, 2) * pow(sin(i), 2))
              / (256 * pow(a, 1.5) * pow(1 - pow(e, 2), 2))
              * (cos(2 * g) * (20 * pow(sin(2 * h), 2) * cos(2 * i) - 30 * cos(4 * h) + 14)
                 + 40 * sin(2 * g) * sin(4 * h) * cos(i))
          + (45 * pow(k, 2) * pow(n3, 2) * pow(a, 3.5) * pow(e, 2))
                / (2048 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                * (4 * cos(2 * g) * pow(sin(i), 2)
                       * (10 * pow(sin(2 * h), 2) * cos(2 * i) - 12 * pow(e, 2) - 15 * cos(4 * h)
                          + 7)
                   + 5
                         * (4 * cos(2 * i) * (-6 * pow(e, 2) + cos(4 * h) + 7) - 40 * pow(e, 2)
                            + 2 * pow(sin(2 * h), 2) * cos(4 * i) - 3 * cos(4 * h) + 35)
                   + 80 * sin(2 * g) * sin(4 * h) * pow(sin(i), 2) * cos(i));
    double H_mp1c
        = (9 * pow(C22_MOON, 2) * pow(GM_MOON, 1.5) * pow(R_MOON, 4) * pow(sin(i), 2) * cos(i))
              / (2 * pow(n3, 2) * pow(a, 13.5) * pow(1 - pow(e, 2), 3.5))
          + (9 * k * C22_MOON * pow(GM_MOON, 0.5) * pow(R_MOON, 2) * pow(sin(i), 2) * cos(i))
                / (16 * pow(a, 1.5) * pow(1 - pow(e, 2), 2))
                * (15 * pow(e, 2) * cos(2 * g) + 6 * pow(e, 2) + 4)
          + (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3.5) * sqrt(1 - pow(e, 2)) * cos(i))
                / (128 * pow(GM_MOON, 0.5))
                * (30 * pow(e, 2) * cos(2 * g) * pow(sin(i), 2) + (17 * pow(e, 2) - 2) * cos(2 * i)
                   + 83 * pow(e, 2) + 2);

    double l_mp1c
        = (45 * k * C22_MOON * pow(R_MOON, 2) * (5 * pow(e, 2) + 2) * pow(sin(i), 2))
              / (256 * pow(a, 2) * pow(1 - pow(e, 2), 2))
              * (sin(2 * g) * (-10 * pow(sin(2 * h), 2) * cos(2 * i) + 15 * cos(4 * h) - 7)
                 + 20 * cos(2 * g) * sin(4 * h) * cos(i))
          - (45 * pow(k, 2) * pow(n3, 2) * pow(a, 3) * sqrt(1 - pow(e, 2)))
                / (4096 * pow(GM_MOON, 0.5))
                * (5
                       * (-8 * sin(4 * h)
                              * ((cos(2 * g) + 7 * cos(4 * g)) * cos(i)
                                 - 2 * pow(sin(g), 2) * (2 * cos(2 * g) + 1) * cos(3 * i))
                          - 16 * pow(sin(g), 3) * cos(g) * pow(sin(2 * h), 2) * cos(4 * i)
                          + sin(4 * g)
                                * (-7 * cos(4 * h) * (4 * cos(2 * i) + 5) - 4 * cos(2 * i) + 3))
                   + 2 * sin(2 * g) * (5 * cos(4 * h) * (4 * cos(2 * i) - 5) - 4 * cos(2 * i) + 9));

    double g_mp1c
        = (9 * pow(C22_MOON, 2) * GM_MOON * pow(R_MOON, 4)
           / (4 * pow(n3, 2) * pow(a, 7) * pow(e * e - 1, 4)) * cos(i) * sin(4 * h))
          - (9 * k * C22_MOON * pow(R_MOON, 2) / (256 * pow(a, 2) * pow(1 - e * e, 5 / 2))
             * (40 * sin(2 * g) * pow(cos(2 * h), 2)
                    * ((3 * pow(e, 2) - 1) * cos(2 * i) - pow(e, 2) + 1)
                + 5 * sin(2 * g) * pow(sin(2 * h), 2)
                      * ((12 - 68 * pow(e, 2)) * cos(2 * i) + (5 - 15 * pow(e, 2)) * cos(4 * i)
                         + 19 * pow(e, 2) - 17)
                + 2 * sin(4 * h) * cos(i)
                      * (25 * (7 * pow(e, 2) - 2) * cos(2 * i) * cos(2 * g)
                         + 25 * (2 - 3 * pow(e, 2)) * cos(2 * g) - 8 * (3 * pow(e, 2) + 2))))
          - (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3) / (8192 * GM_MOON * (1 - pow(e, 2)))
             * (-400 * (3 * pow(e, 2) - 2) * pow(sin(g), 3) * cos(g) * pow(sin(2 * h), 2)
                    * cos(4 * i)
                - 8 * sin(4 * h)
                      * (cos(i)
                             * (36 * pow(e, 4) - 50 * (2 * pow(e, 2) + 1) * cos(2 * g)
                                + 23 * pow(e, 2)
                                - 25 * (4 * pow(e, 4) - 21 * pow(e, 2) + 14) * cos(4 * g) + 16)
                         - 100 * pow(sin(g), 2) * cos(3 * i)
                               * ((3 * pow(e, 2) - 2) * cos(2 * g)))));

    double h_mp1c
        = (-9 * pow(C22_MOON, 2) * GM_MOON * pow(R_MOON, 4) * sin(4 * h) * (cos(2 * i) + 3)
           / (16 * pow(n3, 2) * pow(a, 7) * pow(1 - pow(e, 2), 4)))
          + (9 * k * C22_MOON * pow(R_MOON, 2) / (128 * pow(a, 2) * pow(1 - pow(e, 2), 5.0 / 2.0)))
                * (20 * pow(e, 2) * sin(2 * g) * (5 * cos(4 * h) - 3) * cos(i)
                   + sin(4 * h)
                         * (5 * pow(e, 2) * cos(2 * g) * (7 * cos(2 * i) + 13)
                            - 2 * (3 * pow(e, 2) + 2) * (cos(2 * i) + 3)))
          - (9 * pow(k, 2) * pow(n3, 2) * pow(a, 3) / (1024 * GM_MOON * (1 - pow(e, 2))))
                * (40 * pow(e, 2) * sin(2 * g) * cos(i)
                       * (5 * cos(4 * h) * (pow(e, 2) * cos(2 * g) - 1) - 3 * pow(e, 2) + 3)
                   + sin(4 * h)
                         * (25 * pow(e, 4) * cos(4 * g) * (cos(2 * i) + 3) + 27 * pow(e, 4)
                            - 10 * pow(e, 2) * (3 * pow(e, 2) + 7) * cos(2 * i) * cos(2 * g)
                            - 10 * pow(e, 2) * (13 - 3 * pow(e, 2)) * cos(2 * g) - 14 * pow(e, 2)
                            + (9 * pow(e, 4) + 62 * pow(e, 2) + 4) * cos(2 * i) + 12));

    std::array<double, 6> ret = {l_mp1c, g_mp1c, h_mp1c, L_mp1c, G_mp1c, H_mp1c};
    return ret;
  }

}  // namespace lupnt
