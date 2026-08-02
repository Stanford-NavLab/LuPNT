#include "lupnt/environment/body.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"

namespace lupnt {
  namespace {
    void PopulateBodyFromData(Body& body, const BodyData& data) {
      body.id = data.id;
      body.name = data.name;
      body.GM = data.GM;
      body.R = data.R;
      body.omega = data.omega;
      body.fixed_frame = data.fixed_frame;
      body.inertial_frame = data.inertial_frame;
      body.units = data.units;
    }

    BodyData ScaleBodyData(BodyData data, const UnitSystem& units) {
      data.GM = units.GravitationalParameter(data.GM.val());
      data.R = units.Length(data.R.val());
      data.omega = units.AngularVelocity(data.omega.val());
      data.units = units;
      return data;
    }
  }  // namespace

  /// @brief Create a Body object for the Moon
  /// @return Body object for the Moon
  Body Body::Moon(int n_max, int m_max, std::string gravity_file) {
    return Body::Moon(SI_UNITS, n_max, m_max, gravity_file);
  }

  Body Body::Moon(const UnitSystem& units, int n_max, int m_max, std::string gravity_file) {
    Body moon;
    PopulateBodyFromData(moon, GetBodyData(BodyId::MOON, units));
    moon.use_gravity_field = (n_max > 1);  // zonal J2 (deg-2, m=0) needs only n_max>=2
    if (moon.use_gravity_field)
      moon.gravity_field = ReadHarmonicGravityField<Real>(gravity_file, n_max, m_max, true, units);
    return moon;
  }

  /// @brief Create a Body object for the Earth
  /// @return Body object for the Earth
  Body Body::Earth(int n_max, int m_max, std::string gravity_file) {
    return Body::Earth(SI_UNITS, n_max, m_max, gravity_file);
  }

  Body Body::Earth(const UnitSystem& units, int n_max, int m_max, std::string gravity_file) {
    Body earth;
    PopulateBodyFromData(earth, GetBodyData(BodyId::EARTH, units));
    earth.use_gravity_field = (n_max > 1);  // zonal J2 (deg-2, m=0) needs only n_max>=2
    if (earth.use_gravity_field)
      earth.gravity_field = ReadHarmonicGravityField<Real>(gravity_file, n_max, m_max, true, units);

    return earth;
  }

  /// @brief Create a Body object for the Sun
  /// @return Body object for the Sun
  Body Body::Sun() { return Body::Sun(SI_UNITS); }

  Body Body::Sun(const UnitSystem& units) {
    Body sun;
    PopulateBodyFromData(sun, GetBodyData(BodyId::SUN, units));
    sun.use_gravity_field = false;
    return sun;
  }

  /// @brief Create a Body object for Mars
  /// @return Body object for Mars
  Body Body::Mars(int n_max, int m_max, std::string gravity_file) {
    return Body::Mars(SI_UNITS, n_max, m_max, gravity_file);
  }

  Body Body::Mars(const UnitSystem& units, int n_max, int m_max, std::string gravity_file) {
    Body mars;
    PopulateBodyFromData(mars, GetBodyData(BodyId::MARS, units));
    mars.use_gravity_field = (n_max > 1);  // zonal J2 (deg-2, m=0) needs only n_max>=2
    if (mars.use_gravity_field)
      mars.gravity_field = ReadHarmonicGravityField<Real>(gravity_file, n_max, m_max, true, units);
    return mars;
  }

  /// @brief Create a Body object for Venus
  /// @return Body object for Venus
  Body Body::Venus(int n_max, int m_max, std::string gravity_file) {
    return Body::Venus(SI_UNITS, n_max, m_max, gravity_file);
  }

  Body Body::Venus(const UnitSystem& units, int n_max, int m_max, std::string gravity_file) {
    Body venus;
    PopulateBodyFromData(venus, GetBodyData(BodyId::VENUS, units));
    venus.use_gravity_field = (n_max > 1);  // zonal J2 (deg-2, m=0) needs only n_max>=2
    if (venus.use_gravity_field)
      venus.gravity_field = ReadHarmonicGravityField<Real>(gravity_file, n_max, m_max, true, units);
    return venus;
  }

  Body Body::Jupiter() { return Body::Jupiter(SI_UNITS); }

  Body Body::Jupiter(const UnitSystem& units) {
    Body jupiter;
    PopulateBodyFromData(jupiter, GetBodyData(BodyId::JUPITER, units));
    jupiter.use_gravity_field = false;
    return jupiter;
  }

  Body Body::Saturn() { return Body::Saturn(SI_UNITS); }

  Body Body::Saturn(const UnitSystem& units) {
    Body saturn;
    PopulateBodyFromData(saturn, GetBodyData(BodyId::SATURN, units));
    saturn.use_gravity_field = false;
    return saturn;
  }

  Body Body::Uranus() { return Body::Uranus(SI_UNITS); }

  Body Body::Uranus(const UnitSystem& units) {
    Body uranus;
    PopulateBodyFromData(uranus, GetBodyData(BodyId::URANUS, units));
    uranus.use_gravity_field = false;
    return uranus;
  }

  Body Body::Neptune() { return Body::Neptune(SI_UNITS); }

  Body Body::Neptune(const UnitSystem& units) {
    Body neptune;
    PopulateBodyFromData(neptune, GetBodyData(BodyId::NEPTUNE, units));
    neptune.use_gravity_field = false;
    return neptune;
  }

  Body Body::CreateBody(BodyId body, int n, int m, std::string gravity_file) {
    return Body::CreateBody(body, SI_UNITS, n, m, gravity_file);
  }

  Body Body::CreateBody(BodyId body, const UnitSystem& units, int n, int m,
                        std::string gravity_file) {
    Logger::Debug(fmt::format("Creating body {}", magic_enum::enum_name(body).data()), "Body");
    if (body == BodyId::SUN) return Body::Sun(units);
    if (body == BodyId::EARTH) {
      if (!gravity_file.empty()) return Body::Earth(units, n, m, gravity_file);
      return Body::Earth(units, n, m);
    }
    if (body == BodyId::MOON) {
      if (!gravity_file.empty()) return Body::Moon(units, n, m, gravity_file);
      return Body::Moon(units, n, m);
    }
    if (body == BodyId::MARS) {
      if (!gravity_file.empty()) return Body::Mars(units, n, m, gravity_file);
      return Body::Mars(units, n, m);
    }
    if (body == BodyId::VENUS) {
      if (!gravity_file.empty()) return Body::Venus(units, n, m, gravity_file);
      return Body::Venus(units, n, m);
    }
    if (body == BodyId::JUPITER) return Body::Jupiter(units);
    if (body == BodyId::SATURN) return Body::Saturn(units);
    if (body == BodyId::URANUS) return Body::Uranus(units);
    if (body == BodyId::NEPTUNE) return Body::Neptune(units);
    LUPNT_CHECK(false, "Body not found", "Body");
  }

  BodyData GetBodyData(BodyId id) { return GetBodyData(id, SI_UNITS); }

  BodyData GetBodyData(BodyId id, const UnitSystem& units) {
    BodyData data;
    switch (id) {
      case BodyId::SUN:
        data = {.id = BodyId::SUN,
                .name = "SUN",
                .GM = GM_SUN,
                .R = R_SUN,
                .fixed_frame = Frame::ICRF,
                .inertial_frame = Frame::ICRF,
                .flattening = SUN_F,
                .omega = OMEGA_SUN};
        break;
      case BodyId::MERCURY:
        data = {.id = BodyId::MERCURY,
                .name = "MERCURY",
                .GM = GM_MERCURY,
                .R = R_MERCURY,
                .fixed_frame = Frame::MERCURY_FIXED,
                .inertial_frame = Frame::MERCURY_CI,
                .flattening = MERCURY_F,
                .omega = OMEGA_MERCURY};
        break;
      case BodyId::VENUS:
        data = {.id = BodyId::VENUS,
                .name = "VENUS",
                .GM = GM_VENUS,
                .R = R_VENUS,
                .fixed_frame = Frame::VENUS_FIXED,
                .inertial_frame = Frame::VENUS_CI,
                .flattening = VENUS_F,
                .omega = OMEGA_VENUS};
        break;
      case BodyId::EARTH:
        data = {.id = BodyId::EARTH,
                .name = "EARTH",
                .GM = GM_EARTH,
                .R = R_EARTH,
                .fixed_frame = Frame::ITRF,
                .inertial_frame = Frame::GCRF,
                .flattening = WGS84_F,
                .omega = OMEGA_EARTH};
        break;
      case BodyId::MOON:
        data = {.id = BodyId::MOON,
                .name = "MOON",
                .GM = GM_MOON,
                .R = R_MOON,
                .fixed_frame = Frame::MOON_PA,
                .inertial_frame = Frame::MOON_CI,
                .flattening = MOON_F,
                .omega = OMEGA_MOON};
        break;
      case BodyId::MARS:
        data = {.id = BodyId::MARS,
                .name = "MARS",
                .GM = GM_MARS,
                .R = R_MARS,
                .fixed_frame = Frame::MARS_FIXED,
                .inertial_frame = Frame::MARS_CI,
                .flattening = MARS_F,
                .omega = OMEGA_MARS};
        break;
      case BodyId::JUPITER:
        data = {.id = BodyId::JUPITER,
                .name = "JUPITER",
                .GM = GM_JUPITER,
                .R = R_JUPITER,
                .fixed_frame = Frame::JUPITER_FIXED,
                .inertial_frame = Frame::JUPITER_CI,
                .flattening = JUPITER_F,
                .omega = OMEGA_JUPITER};
        break;
      case BodyId::SATURN:
        data = {.id = BodyId::SATURN,
                .name = "SATURN",
                .GM = GM_SATURN,
                .R = R_SATURN,
                .fixed_frame = Frame::SATURN_FIXED,
                .inertial_frame = Frame::SATURN_CI,
                .flattening = SATURN_F,
                .omega = OMEGA_SATURN};
        break;
      case BodyId::URANUS:
        data = {.id = BodyId::URANUS,
                .name = "URANUS",
                .GM = GM_URANUS,
                .R = R_URANUS,
                .fixed_frame = Frame::URANUS_FIXED,
                .inertial_frame = Frame::URANUS_CI,
                .flattening = URANUS_F,
                .omega = OMEGA_URANUS};
        break;
      case BodyId::NEPTUNE:
        data = {.id = BodyId::NEPTUNE,
                .name = "NEPTUNE",
                .GM = GM_NEPTUNE,
                .R = R_NEPTUNE,
                .fixed_frame = Frame::NEPTUNE_FIXED,
                .inertial_frame = Frame::NEPTUNE_CI,
                .flattening = NEPUTUNE_F,
                .omega = OMEGA_NEPTUNE};
        break;
      default: throw std::runtime_error("Body not found");
    }
    return ScaleBodyData(data, units);
  }

  double GetBodyRadius(BodyId body) { return GetBodyRadius(body, SI_UNITS); }

  double GetBodyRadius(BodyId body, const UnitSystem& units) {
    BodyData data = GetBodyData(body, units);
    return data.R.val();
  }

  double GetBodyGM(BodyId body) { return GetBodyGM(body, SI_UNITS); }

  double GetBodyGM(BodyId body, const UnitSystem& units) {
    BodyData data = GetBodyData(body, units);
    return data.GM.val();
  }

  double GetBodyFlattening(BodyId body) {
    BodyData data = GetBodyData(body);
    return data.flattening.val();
  }

  double GetBodyOmega(BodyId body) { return GetBodyOmega(body, SI_UNITS); }

  double GetBodyOmega(BodyId body, const UnitSystem& units) {
    BodyData data = GetBodyData(body, units);
    return data.omega.val();
  }

  std::string GetBodyName(BodyId body) {
    BodyData data = GetBodyData(body);
    return data.name;
  }

  Frame GetInertialFrameName(BodyId body) {
    BodyData data = GetBodyData(body);
    return data.inertial_frame;
  }

  Frame GetBodyFixedFrameName(BodyId body) {
    BodyData data = GetBodyData(body);
    return data.fixed_frame;
  }

  Body CreateBody(BodyId body, int n, int m) { return CreateBody(body, SI_UNITS, n, m); }

  Body CreateBody(BodyId body, const UnitSystem& units, int n, int m) {
    return Body::CreateBody(body, units, n, m);
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
    return ReadHarmonicGravityField<T>(filename, n, m, normalized, SI_UNITS);
  }

  template <typename T> GravityField<T> ReadHarmonicGravityField(const std::string& filename, int n,
                                                                 int m, bool normalized,
                                                                 const UnitSystem& units) {
    GravityField<T> gravity_field;
    gravity_field.units = units;
    std::filesystem::path filepath = GetFilePath(filename);
    std::ifstream file = OpenFile<std::ifstream>(filepath);

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
        gravity_field.GM = units.GravitationalParameter(GM);
        gravity_field.R = units.Length(r);
        break;
      }
    }
    gravity_field.n = n;
    gravity_field.m = m;

    // Initialize Eigen matrices with the specified maxN and maxM
    if (n > gravity_field.n_max || m > gravity_field.m_max)
      throw std::runtime_error("Degree and order exceed maximum values");

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

  // Explicit instantiations so the public ReadHarmonicGravityField overloads
  // link from other translation units (e.g. the Orekit force-model tests),
  // not only from the Body::Earth/Moon calls inside this file.
  template GravityField<Real> ReadHarmonicGravityField<Real>(const std::string&, int, int, bool);
  template GravityField<Real> ReadHarmonicGravityField<Real>(const std::string&, int, int, bool,
                                                             const UnitSystem&);
  template GravityField<double> ReadHarmonicGravityField<double>(const std::string&, int, int,
                                                                 bool);
  template GravityField<double> ReadHarmonicGravityField<double>(const std::string&, int, int, bool,
                                                                 const UnitSystem&);

}  // namespace lupnt
