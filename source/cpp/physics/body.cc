#include "lupnt/physics/body.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "lupnt/core/user_file_path.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  /// @brief Create a Body object for the Moon
  /// @return Body object for the Moon
  Body Body::Moon(int n_max, int m_max, std::string gravity_file) {
    Body moon;
    moon.name = "MOON";
    moon.id = NaifId::MOON;
    moon.fixed_frame = Frame::MOON_PA;
    moon.inertial_frame = Frame::MOON_CI;
    moon.GM = GM_MOON;
    moon.R = R_MOON;
    moon.gravity_field = ReadHarmonicGravityField(gravity_file, n_max, m_max, true);
    return moon;
  }

  /// @brief Create a Body object for the Earth
  /// @return Body object for the Earth
  Body Body::Earth(int n_max, int m_max, std::string gravity_file) {
    Body earth;
    earth.name = "EARTH";
    earth.id = NaifId::EARTH;
    earth.fixed_frame = Frame::ITRF;
    earth.inertial_frame = Frame::GCRF;
    earth.GM = GM_EARTH;
    earth.R = R_EARTH;
    earth.gravity_field = ReadHarmonicGravityField(gravity_file, n_max, m_max, true);

    return earth;
  }

  /// @brief Create a Body object for the Sun
  /// @return Body object for the Sun
  Body Body::Sun() {
    Body sun;
    sun.name = "SUN";
    sun.id = NaifId::SUN;
    sun.inertial_frame = Frame::ICRF;
    sun.GM = GM_SUN;
    sun.R = 696342.0;
    sun.use_gravity_field = false;
    return sun;
  }

  /// @brief Create a Body object for Mars
  /// @return Body object for Mars
  Body Body::Mars(int n_max, int m_max, std::string gravity_file) {
    Body mars;
    mars.name = "MARS";
    mars.id = NaifId::MARS;
    mars.fixed_frame = Frame::MARS_FIXED;
    mars.GM = 0.4282837566395650E+05;
    mars.R = 0.3396000000000000E+04;
    mars.gravity_field = ReadHarmonicGravityField(gravity_file, n_max, m_max, true);
    return mars;
  }

  /// @brief Create a Body object for Venus
  /// @return Body object for Venus
  Body Body::Venus(int n_max, int m_max, std::string gravity_file) {
    Body venus;
    venus.name = "VENUS";
    venus.id = NaifId::VENUS;
    venus.fixed_frame = Frame::VENUS_FIXED;
    venus.GM = 0.3248585920790000E+06;
    venus.R = 0.6051000000000000E+04;
    venus.gravity_field = ReadHarmonicGravityField(gravity_file, n_max, m_max, true);
    return venus;
  }

  /// @brief Kronecker delta function
  /// @param i
  /// @param j
  /// @return
  double kron(int i, int j) { return (i == j) ? 1 : 0; }

  /// @brief Compute the factorial product (n-m)!/(n+m)!
  /// @param n
  /// @param m
  /// @return
  double factprod(int n, int m) {
    double f = 1.0;
    for (int i = n - m + 1; i <= n + m; i++) {
      f /= i;
    }
    return f;
  }

  /// @brief Read a harmonic gravity field from a file
  /// @param filename Harmonic gravity field filename
  /// @param n Degree of the spherical harmonics expansion
  /// @param m Order of the spherical harmonics expansion
  /// @param normalized Whether the coefficients are normalized
  /// @return Gravity field object
  GravityField ReadHarmonicGravityField(const std::string& filename, int n_max, int m_max,
                                        bool normalized) {
    GravityField gravity_field;
    std::filesystem::path filepath = GetFilePath(filename);
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    // Read header lines
    std::string line;
    while (std::getline(file, line)) {
      if (line.find("POTFIELD") != std::string::npos) {
        std::string potfield = line.substr(0, 8);
        int n_max = std::stoi(line.substr(8, 3));
        int m_max = std::stoi(line.substr(11, 3));

        std::istringstream iss(line.substr(14));
        double GM, r, dummyFactor;
        iss >> GM >> r >> dummyFactor;
        gravity_field.n_max = n_max;
        gravity_field.m_max = m_max;
        gravity_field.GM = GM * pow(KM_M, 3);
        gravity_field.R = r * KM_M;
        break;
      }
    }

    // Initialize Eigen matrices with the specified maxN and maxM
    assert(n_max <= gravity_field.n_max && m_max <= gravity_field.m_max);
    gravity_field.CS = Eigen::MatrixXd::Zero(n_max + 1, m_max + 1);
    gravity_field.CS(0, 0) = 1.0;  // C00 = 1.0
    // Read coefficient lines
    while (std::getline(file, line)) {
      if (line.find("RECOEF") != std::string::npos) {
        std::string recoef = line.substr(0, 8);
        int n = std::stoi(line.substr(8, 3));
        int m = std::stoi(line.substr(11, 3));
        std::istringstream iss(line.substr(14));
        double cnm, snm = 0.0;
        iss >> cnm;
        if (n >= n_max + 1) break;
        if (m >= m_max + 1) continue;

        if (m == 0) {
          double N = (normalized) ? sqrt(2 * n + 1) : 1.0;
          gravity_field.CS(n, m) = N * cnm;
        } else {
          double N = (normalized) ? sqrt((2 - kron(0, m)) * (2 * n + 1) * factprod(n, m)) : 1.0;
          iss >> snm;
          gravity_field.CS(n, m) = N * cnm;
          gravity_field.CS(m - 1, n) = N * snm;
        }
      } else if (line.find("END") != std::string::npos) {
        break;
      }
    }

    file.close();
    return gravity_field;
  }

}  // namespace lupnt
