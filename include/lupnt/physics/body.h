/**
 * @file agent.h
 * @author Stanford NAV LAB
 * @brief List of agents
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string>

#include "lupnt/core/constants.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  struct GravityField {
    int n_max, m_max;  // Maximum degree and order
    int n, m;          // Degree and order

    Real GM;  // Gravitational constant [km^3/s^2]
    Real R;   // Reference radius [km]
    MatX CS;  // Unnormalized coefficients
  };

  struct Body {
    NaifId id;
    std::string name;

    Real GM;
    Real R;
    int n, m;

    Frame fixed_frame;
    Frame inertial_frame;

    bool use_gravity_field;
    GravityField gravity_field;

    static Body Sun();
    static Body Earth(int n_max = 0, int m_max = 0, std::string gravity_file = "EGM96.cof");
    static Body Moon(int n_max = 0, int m_max = 0, std::string gravity_file = "grgm900c.cof");
    static Body Venus(int n_max = 0, int m_max = 0, std::string gravity_file = "MGN75HSAAP.cof");
    static Body Mars(int n_max = 0, int m_max = 0, std::string gravity_file = "GMM1.cof");
  };

  GravityField ReadHarmonicGravityField(const std::string& filename, int n, int m, bool normalized);
}  // namespace lupnt
