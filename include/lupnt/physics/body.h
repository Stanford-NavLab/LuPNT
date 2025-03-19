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

  struct BodyData {
    NaifId id;
    std::string name;
    Real GM;
    Real R;
    Frame fixed_frame;
    Frame inertial_frame;
  };

  template <typename T = double> struct GravityField {
    int n_max, m_max;  // Maximum degree and order
    int n, m;          // Degree and order

    T GM;                            // Gravitational constant [km^3/s^2]
    T R;                             // Reference radius [km]
    Matrix<T, Dynamic, Dynamic> CS;  // Unnormalized coefficients
  };

  template <typename T = double> struct BodyT {
    NaifId id;
    std::string name;

    T GM;
    T R;
    int n, m;

    Frame fixed_frame;
    Frame inertial_frame;

    bool use_gravity_field;
    GravityField<T> gravity_field;

    static BodyT Sun();
    static BodyT Earth(int n = 0, int m = 0, std::string gravity_file = "EGM96.cof");
    static BodyT Moon(int n = 0, int m = 0, std::string gravity_file = "grgm900c.cof");
    static BodyT Venus(int n = 0, int m = 0, std::string gravity_file = "MGN75HSAAP.cof");
    static BodyT Mars(int n = 0, int m = 0, std::string gravity_file = "GMM1.cof");
  };

  using Body = BodyT<double>;

  template <typename T> GravityField<T> ReadHarmonicGravityField(const std::string& filename, int n,
                                                                 int m, bool normalized);

  BodyData GetBodyData(NaifId id);
  double GetBodyRadius(NaifId body);
  std::string GetBodyName(NaifId body);
  Frame GetInertialFrameName(NaifId body);
  Frame GetBodyFixedFrameName(NaifId body);

}  // namespace lupnt
