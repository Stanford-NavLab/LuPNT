#include "lupnt/physics/body.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "lupnt/physics/frame_converter.h"

namespace lupnt {
Body Body::Moon(int n_max, int m_max) {
  Body moon;
  moon.name = "MOON";
  moon.id = NaifId::MOON;
  moon.spherical_harmonics = n_max > 0 || m_max > 0;
  moon.normalized = true;
  moon.n_max = n_max;
  moon.m_max = m_max;
  moon.fixed_frame = Frame::MOON_PA;

  BodyData bd = GetBodyData(moon.id);
  moon.GM = bd.GM;
  moon.R = bd.R;

  if (moon.spherical_harmonics)
    std::tie(moon.Cnm, moon.Snm) = LoadGravityCoefficients(bd, n_max);
  return moon;
}

Body Body::Earth(int n_max, int m_max) {
  Body earth;
  earth.name = "EARTH";
  earth.id = NaifId::EARTH;
  earth.spherical_harmonics = n_max > 0 || m_max > 0;
  earth.normalized = true;
  earth.n_max = n_max;
  earth.m_max = m_max;
  earth.fixed_frame = Frame::ITRF;

  BodyData bd = GetBodyData(earth.id);
  earth.GM = bd.GM;
  earth.R = bd.R;

  if (earth.spherical_harmonics)
    std::tie(earth.Cnm, earth.Snm) = LoadGravityCoefficients(bd, n_max);
  return earth;
}

Body Body::Sun() {
  Body sun;
  sun.name = "SUN";
  sun.id = NaifId::SUN;
  sun.spherical_harmonics = false;
  sun.normalized = false;
  sun.fixed_frame = Frame::ICRF;

  BodyData bd = GetBodyData(sun.id);
  sun.GM = bd.GM;
  sun.R = bd.R;
  return sun;
}

Body Body::Mars(int n_max, int m_max) {
  Body mars;
  mars.name = "MARS";
  mars.id = NaifId::MARS;
  mars.spherical_harmonics = n_max > 0 || m_max > 0;
  mars.normalized = true;
  mars.n_max = n_max;
  mars.m_max = m_max;
  mars.fixed_frame = Frame::MARS_FIXED;

  BodyData bd = GetBodyData(mars.id);
  mars.GM = bd.GM;
  mars.R = bd.R;

  if (mars.spherical_harmonics)
    std::tie(mars.Cnm, mars.Snm) = LoadGravityCoefficients(bd, n_max);
  return mars;
}

Body Body::Venus(int n_max, int m_max) {
  Body venus;
  venus.name = "VENUS";
  venus.id = NaifId::VENUS;
  venus.spherical_harmonics = n_max > 0 || m_max > 0;
  venus.normalized = true;
  venus.n_max = n_max;
  venus.m_max = m_max;
  venus.fixed_frame = Frame::VENUS_FIXED;

  BodyData bd = GetBodyData(venus.id);
  venus.GM = bd.GM;
  venus.R = bd.R;

  if (venus.spherical_harmonics)
    std::tie(venus.Cnm, venus.Snm) = LoadGravityCoefficients(bd, n_max);
  return venus;
}

BodyData GetBodyData(const NaifId bodyID) {
  BodyData bd;

  switch (bodyID) {
    case NaifId::VENUS:  // VENUS (180x180)
      bd.filepath = "shgj180u.a01";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 0.3248585920790000E+06;
      bd.R = 0.6051000000000000E+04;
      break;

    case NaifId::EARTH:  // EARTH (50x50)
      bd.filepath = "GGM02C.GEO";
      bd.headerlines = 3;
      bd.delimiter = "";
      bd.GM = 398600.44150;
      bd.R = 6378136.30 * 1e-3;
      break;

    case NaifId::MOON:  // MOON (660x660)
      bd.filepath = "gggrx_0660pm_sha.tab";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 4.9028001224453001e+03;
      bd.R = 1.7380000000000000e+03;
      break;

    case NaifId::MARS:  // MARS (120x120)
      bd.filepath = "jgmro_120f_sha.tab";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 0.4282837566395650E+05;
      bd.R = 0.3396000000000000E+04;
      break;

    case NaifId::SUN:  // SUN (10x10)
      bd.filepath = "shgm405c.bsp";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 132712440041.9394;
      bd.R = 696342.0;
      break;

    default:
      throw std::runtime_error(
          "Spherical Harmonics coefficients not available!");
  }

  return bd;
}

std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

void ReadData(const std::string &filepath, int N, int headerlines,
              const std::string &delimiter, std::vector<int> &idN,
              std::vector<int> &idM, std::vector<double> &C,
              std::vector<double> &S, std::vector<double> &sigC,
              std::vector<double> &sigS) {
  std::ifstream file(filepath);

  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filepath);
  }

  // Skip header lines
  std::string line;
  for (int i = 0; i < headerlines; ++i) {
    std::getline(file, line);
  }

  // Read and parse data lines
  for (int i = 0; i < N && std::getline(file, line); ++i) {
    std::vector<std::string> tokens = split(line, ',');

    idN[i] = std::stoi(tokens[0]);
    idM[i] = std::stoi(tokens[1]);
    C[i] = std::stod(tokens[2]);
    S[i] = std::stod(tokens[3]);
    sigC[i] = std::stod(tokens[4]);
    sigS[i] = std::stod(tokens[5]);
  }

  file.close();
}

std::tuple<VecXd, VecXd> LoadGravityCoefficients(BodyData bd, int nmax) {
  // read text file
  int N = nmax * nmax + 10;

  std::vector<int> idN(N), idM(N);
  std::vector<double> C(N), S(N), sigC(N), sigS(N);

  std::filesystem::path harmonicsPath = GetDataPath() / "spherical_harmonics";
  ReadData(harmonicsPath / bd.filepath, N, bd.headerlines, bd.delimiter, idN,
           idM, C, S, sigC, sigS);

  VecXd Cnm = VecXd::Zero(nmax + 2, nmax + 2);
  VecXd Snm = VecXd::Zero(nmax + 2, nmax + 2);

  Cnm(0, 0) = 1.0;
  Cnm(1, 0) = 0.0;
  Cnm(1, 1) = 0.0;

  Snm(0, 0) = 0.0;
  Snm(1, 0) = 0.0;
  Snm(1, 1) = 0.0;

  for (int k = 0; k < (((nmax + 1) * (nmax + 2) / 2) - 3); ++k) {
    Cnm(idN[k], idM[k]) = C[k];
    Snm(idN[k], idM[k]) = S[k];
  }

  return {Cnm, Snm};
}

}  // namespace lupnt