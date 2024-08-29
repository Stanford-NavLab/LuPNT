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
  template struct BodyT<double>;
  template struct BodyT<Real>;

  /// @brief Create a Body object for the Moon
  /// @return Body object for the Moon
  template <typename T> BodyT<T> BodyT<T>::Moon(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> moon;
    moon.name = "MOON";
    moon.id = NaifId::MOON;
    moon.fixed_frame = Frame::MOON_PA;
    moon.inertial_frame = Frame::MOON_CI;
    moon.GM = GM_MOON;
    moon.R = R_MOON;
    moon.use_gravity_field = (n_max > 1 && m_max > 1);
    if (moon.use_gravity_field)
      moon.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);
    return moon;
  }
  template BodyT<double> BodyT<double>::Moon(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Moon(int n_max, int m_max, std::string gravity_file);

  /// @brief Create a Body object for the Earth
  /// @return Body object for the Earth
  template <typename T> BodyT<T> BodyT<T>::Earth(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> earth;
    earth.name = "EARTH";
    earth.id = NaifId::EARTH;
    earth.fixed_frame = Frame::ITRF;
    earth.inertial_frame = Frame::GCRF;
    earth.GM = GM_EARTH;
    earth.R = R_EARTH;
    earth.use_gravity_field = (n_max > 1 && m_max > 1);
    if (earth.use_gravity_field)
      earth.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);

    return earth;
  }
  template BodyT<double> BodyT<double>::Earth(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Earth(int n_max, int m_max, std::string gravity_file);

  /// @brief Create a Body object for the Sun
  /// @return Body object for the Sun
  template <typename T> BodyT<T> BodyT<T>::Sun() {
    BodyT<T> sun;
    sun.name = "SUN";
    sun.id = NaifId::SUN;
    sun.inertial_frame = Frame::ICRF;
    sun.GM = GM_SUN;
    sun.R = R_SUN;
    sun.use_gravity_field = false;
    return sun;
  }
  template BodyT<double> BodyT<double>::Sun();
  template BodyT<Real> BodyT<Real>::Sun();

  /// @brief Create a Body object for Mars
  /// @return Body object for Mars
  template <typename T> BodyT<T> BodyT<T>::Mars(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> mars;
    mars.name = "MARS";
    mars.id = NaifId::MARS;
    mars.fixed_frame = Frame::MARS_FIXED;
    mars.GM = GM_MARS;
    mars.R = R_MARS;
    mars.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);
    return mars;
  }
  template BodyT<double> BodyT<double>::Mars(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Mars(int n_max, int m_max, std::string gravity_file);

  /// @brief Create a Body object for Venus
  /// @return Body object for Venus
  template <typename T> BodyT<T> BodyT<T>::Venus(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> venus;
    venus.name = "VENUS";
    venus.id = NaifId::VENUS;
    venus.fixed_frame = Frame::VENUS_FIXED;
    venus.GM = GM_MARS;
    venus.R = R_VENUS;
    venus.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);
    return venus;
  }
  template BodyT<double> BodyT<double>::Venus(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Venus(int n_max, int m_max, std::string gravity_file);

  BodyData GetBodyData(NaifId id) {
    switch (id) {
      case NaifId::SUN: return {NaifId::SUN, "SUN", GM_SUN, 696342.0, Frame::ICRF, Frame::ICRF};
      case NaifId::MERCURY:
        return {NaifId::MERCURY,      "MERCURY",        GM_MERCURY, R_MERCURY,
                Frame::MERCURY_FIXED, Frame::MERCURY_CI};
      case NaifId::VENUS:
        return {NaifId::VENUS, "VENUS", GM_VENUS, R_VENUS, Frame::VENUS_FIXED, Frame::VENUS_CI};
      case NaifId::EARTH:
        return {NaifId::EARTH, "EARTH", GM_EARTH, R_EARTH, Frame::ITRF, Frame::GCRF};
      case NaifId::MOON:
        return {NaifId::MOON, "MOON", GM_MOON, R_MOON, Frame::MOON_PA, Frame::MOON_CI};
      case NaifId::MARS:
        return {NaifId::MARS, "MARS", GM_MARS, R_MARS, Frame::MARS_FIXED, Frame::MARS_CI};
      case NaifId::JUPITER:
        return {NaifId::JUPITER,      "JUPITER",        GM_JUPITER, R_JUPITER,
                Frame::JUPITER_FIXED, Frame::JUPITER_CI};
      case NaifId::SATURN:
        return {NaifId::SATURN,      "SATURN",        GM_SATURN, R_SATURN,
                Frame::SATURN_FIXED, Frame::SATURN_CI};
      case NaifId::URANUS:
        return {NaifId::URANUS,      "URANUS",        GM_URANUS, R_URANUS,
                Frame::URANUS_FIXED, Frame::URANUS_CI};
      case NaifId::NEPTUNE:
        return {NaifId::NEPTUNE,      "NEPTUNE",        GM_NEPTUNE, R_NEPTUNE,
                Frame::NEPTUNE_FIXED, Frame::NEPTUNE_CI};
      default: break;
    }
    throw std::runtime_error("Body not found");
    return {NaifId::SUN, "SUN", GM_SUN, 696342.0, Frame::ICRF, Frame::ICRF};
  }

  double GetBodyRadius(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.R.val();
  }

  std::string GetBodyName(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.name;
  }

  Frame GetInertialFrameName(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.inertial_frame;
  }

  Frame GetBodyFixedFrameName(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.fixed_frame;
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
  template <typename T> GravityField<T> ReadHarmonicGravityField(const std::string& filename, int n,
                                                                 int m, bool normalized) {
    GravityField<T> gravity_field;
    std::filesystem::path filepath = GetFilePath(filename);
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    // Read header lines
    std::string line;
    while (std::getline(file, line)) {
      if (line.find("POTFIELD") != std::string::npos) {
        std::string potfield = line.substr(0, 8);
        int n_max_in = std::stoi(line.substr(8, 3));
        int m_max_in = std::stoi(line.substr(11, 3));

        std::istringstream iss(line.substr(14));
        double dummy1, GM, r, dummy2;
        iss >> dummy1 >> GM >> r >> dummy2;
        gravity_field.n_max = n_max_in;
        gravity_field.m_max = m_max_in;
        gravity_field.GM = GM * pow(KM_M, 3);
        gravity_field.R = r * KM_M;
        break;
      }
    }
    gravity_field.n = n;
    gravity_field.m = m;

    // Initialize Eigen matrices with the specified maxN and maxM
    assert(n <= gravity_field.n_max && m <= gravity_field.m_max);
    gravity_field.CS = Eigen::MatrixXd::Zero(n + 1, m + 1);
    gravity_field.CS(0, 0) = 1.0;  // C00 = 1.0
    // Read coefficient lines
    while (std::getline(file, line)) {
      if (line.find("RECOEF") != std::string::npos) {
        std::string recoef = line.substr(0, 8);
        int n_in = std::stoi(line.substr(8, 3));
        int m_in = std::stoi(line.substr(11, 3));
        std::istringstream iss(line.substr(14));
        double cnm, snm = 0.0;
        iss >> cnm;
        if (n_in >= n + 1) break;
        if (m_in >= m + 1) continue;

        if (m_in == 0) {
          double N = (normalized) ? sqrt(2 * n_in + 1) : 1.0;
          gravity_field.CS(n_in, m_in) = N * cnm;
        } else {
          double N = (normalized)
                         ? sqrt((2 - kron(0, m_in)) * (2 * n_in + 1) * factprod(n_in, m_in))
                         : 1.0;
          iss >> snm;
          double C = N * cnm;
          double S = N * snm;
          gravity_field.CS(n_in, m_in) = C;
          gravity_field.CS(m_in - 1, n_in) = S;
        }
      } else if (line.find("END") != std::string::npos) {
        break;
      }
    }

    file.close();
    return gravity_field;
  }

}  // namespace lupnt
