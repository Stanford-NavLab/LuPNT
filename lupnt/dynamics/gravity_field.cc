/**
 * @file GravityField.cpp
 * @author Stanford NAV LAB
 * @brief Gravity Field Implementations
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "gravity_field.h"

namespace lupnt {

std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

BodyData GetBodyData(const BodyId bodyID) {
  BodyData bd;

  switch (bodyID) {
    case BodyId::VENUS:  // VENUS (180x180)
      bd.filepath = "shgj180u.a01";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 0.3248585920790000E+06;
      bd.R = 0.6051000000000000E+04;
      break;

    case BodyId::EARTH:  // EARTH (50x50)
      bd.filepath = "GGM02C.GEO";
      bd.headerlines = 3;
      bd.delimiter = "";
      bd.GM = 398600.44150;
      bd.R = 6378136.30 * 1e-3;
      break;

    case BodyId::MOON:  // MOON (660x660)
      bd.filepath = "gggrx_0660pm_sha.tab";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 4.9028001224453001e+03;
      bd.R = 1.7380000000000000e+03;
      break;

    case BodyId::MARS:  // MARS (120x120)
      bd.filepath = "jgmro_120f_sha.tab";
      bd.headerlines = 3;
      bd.delimiter = ",";
      bd.GM = 0.4282837566395650E+05;
      bd.R = 0.3396000000000000E+04;
      break;

    default:
      throw std::runtime_error(
          "Spherical Harmonics coefficients not available!");
  }

  return bd;
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

std::tuple<MatrixXd, MatrixXd> LoadGravityCoefficients(BodyData bd, int nmax) {
  // read text file
  int N = nmax * nmax + 10;

  std::vector<int> idN(N), idM(N);
  std::vector<double> C(N), S(N), sigC(N), sigS(N);

  std::filesystem::path harmonicsPath = GetDataPath() / "spherical_harmonics";
  ReadData(harmonicsPath / bd.filepath, N, bd.headerlines, bd.delimiter, idN,
           idM, C, S, sigC, sigS);

  MatrixXd Cnm = MatrixXd::Zero(nmax + 2, nmax + 2);
  MatrixXd Snm = MatrixXd::Zero(nmax + 2, nmax + 2);

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

std::tuple<double, double> spharm_vwmm(int m_in, double Vm_1m_1, double Wm_1m_1,
                                       const Vector3d &x_R, double Re) {
  double m = m_in;

  Vector3d a = x_R;
  int kdel0, kdel1;
  if (m_in > 1) {
    kdel0 = 2;
    kdel1 = 2;
  } else if (m_in == 1) {
    kdel0 = 2;
    kdel1 = 1;
  } else {
    kdel0 = 1;
    kdel1 = 2;
  }

  double r2 = a.dot(a);
  double Vmm = sqrt(kdel0 / kdel1 * (2 * m + 1) / (2 * m)) *
               (a(0) * Re / r2 * Vm_1m_1 - a(1) * Re / r2 * Wm_1m_1);
  double Wmm = sqrt(kdel0 / kdel1 * (2 * m + 1) / (2 * m)) *
               (a(0) * Re / r2 * Wm_1m_1 + a(1) * Re / r2 * Vm_1m_1);
  return {Vmm, Wmm};
}

std::tuple<double, double> spharm_vwm1m(int m_in, double Vmm, double Wmm,
                                        const Vector3d &x_R, double Re) {
  double m = m_in;

  double r2 = x_R.dot(x_R);
  double Vm1m = sqrt(2 * m + 3) * x_R(2) * Re / r2 * Vmm;
  double Wm1m = sqrt(2 * m + 3) * x_R(2) * Re / r2 * Wmm;
  return {Vm1m, Wm1m};
}

std::tuple<double, double> spharm_vwnm(int n_in, int m_in, double Vn_1m,
                                       double Vn_2m, double Wn_1m, double Wn_2m,
                                       const Vector3d &x_R, double Re) {
  double n = n_in;
  double m = m_in;

  Vector3d a = x_R;
  double r2 = a.dot(a);
  double anm = sqrt(((2 * n + 1) * (2 * n - 1)) / ((n + m) * (n - m)));
  double bnm = -sqrt(((2 * n + 1) * (n + m - 1) * (n - m - 1)) /
                     ((2 * n - 3) * (n + m) * (n - m)));

  double Vnm = anm * a(2) * Re / r2 * Vn_1m + bnm * (Re * Re) / r2 * Vn_2m;
  double Wnm = anm * a(2) * Re / r2 * Wn_1m + bnm * (Re * Re) / r2 * Wn_2m;
  return {Vnm, Wnm};
}

Vector3d Facc_j(const Vector3d &facc_R, const Matrix3d &Ur2j) {
  Vector3d acc_j = Ur2j * facc_R;
  return acc_j;
}

std::tuple<Vector3d, Vector3d> spharm_dvwdx(int n_in, int m_in, double Vn1m,
                                            double Vn1m1, double Vn1m_1,
                                            double Wn1m, double Wn1m1,
                                            double Wn1m_1, double Re) {
  Vector3d dVdX = Vector3d::Zero();
  Vector3d dWdX = Vector3d::Zero();

  double Ca;
  int kdel;

  double n = n_in;
  double m = m_in;

  if (m_in == 0) {
    Ca =
        -1.0 / Re * sqrt(((2 * n + 1) * (n + 2) * (n + 1)) / (2 * (2 * n + 3)));

    dVdX(0) = Ca * Vn1m1;
    dVdX(1) = Ca * Wn1m1;
  } else {
    if (m_in != 1)
      kdel = 2;
    else
      kdel = 1;

    Ca = 0.5 / Re * sqrt((2 * n + 1) / (2 * n + 3));
    double Cb = sqrt((n + m + 2) * (n + m + 1));
    double Cc = sqrt(2 * (n - m + 2) * (n - m + 1) / kdel);

    dVdX(0) = -Ca * (Cb * Vn1m1 - Cc * Vn1m_1);
    dVdX(1) = -Ca * (Cb * Wn1m1 + Cc * Wn1m_1);
    dWdX(0) = -Ca * (Cb * Wn1m1 - Cc * Wn1m_1);
    dWdX(1) = Ca * (Cb * Vn1m1 + Cc * Vn1m_1);
  }

  double Cd =
      -1.0 / Re * sqrt((2 * n + 1) * (n + m + 1) * (n - m + 1) / (2 * n + 3));
  dVdX(2) = Cd * Vn1m;
  dWdX(2) = Cd * Wn1m;

  return {dVdX, dWdX};
}

std::tuple<Matrix3d, Matrix3d> spharm_d2vwdx2(int n_in, int m_in, double Vn2m,
                                              double Vn2m1, double Vn2m2,
                                              double Vn2m_1, double Vn2m_2,
                                              double Wn2m, double Wn2m1,
                                              double Wn2m2, double Wn2m_1,
                                              double Wn2m_2, double Re) {
  double n = n_in;
  double m = m_in;

  Matrix3d d2VdX2 = Matrix3d::Zero();
  Matrix3d d2WdX2 = Matrix3d::Zero();
  VectorXd cf(11);
  Vector2d kdel;

  cf(0) = 1 / (Re * Re) * sqrt((2 * n + 1) / (2 * n + 5));
  cf(1) = sqrt((n + m + 2) * (n + m + 1) * (n - m + 2) * (n - m + 1));

  d2VdX2(2, 2) = cf(0) * cf(1) * Vn2m;
  d2WdX2(2, 2) = cf(0) * cf(1) * Wn2m;

  if (m_in == 0) {
    cf(2) = sqrt((n + 4) * (n + 3) * (n + 2) * (n + 1) / 2);
    cf(3) = (n + 2) * (n + 1);
    cf(4) = sqrt((n + 1) * (n + 3) * (n + 2) * (n + 1) / 2);

    d2VdX2(0, 0) = cf(0) / 2 * (cf(2) * Vn2m2 - cf(3) * Vn2m);
    d2VdX2(0, 1) = cf(0) / 2 * (cf(2) * Wn2m2);
    d2VdX2(0, 2) = cf(0) * cf(4) * Vn2m1;
    d2VdX2(1, 1) = cf(0) / 2 * (-cf(2) * Vn2m2 - cf(3) * Vn2m);
    d2VdX2(1, 2) = cf(0) * cf(4) * Wn2m1;

    d2WdX2(0, 0) = 0;
    d2WdX2(0, 1) = 0;
    d2WdX2(0, 2) = 0;
    d2WdX2(1, 1) = 0;
    d2WdX2(1, 2) = 0;
  } else {
    if (m_in != 1)
      kdel(0) = 2;
    else
      kdel(0) = 1;

    cf(5) = sqrt((n - m + 1) * (n + m + 3) * (n + m + 2) * (n + m + 1));
    cf(6) = sqrt(2 / kdel(0) * (n + m + 1) * (n - m + 3) * (n - m + 2) *
                 (n - m + 1));

    d2VdX2(0, 2) = cf(0) / 2 * (cf(5) * Vn2m1 - cf(6) * Vn2m_1);
    d2VdX2(1, 2) = cf(0) / 2 * (cf(5) * Wn2m1 + cf(6) * Wn2m_1);

    d2WdX2(0, 2) = cf(0) / 2 * (cf(5) * Wn2m1 - cf(6) * Wn2m_1);
    d2WdX2(1, 2) = cf(0) / 2 * (-cf(5) * Vn2m1 - cf(6) * Vn2m_1);

    if (m_in == 1) {
      cf(7) = sqrt((n + 5) * (n + 4) * (n + 3) * (n + 2));
      cf(8) = sqrt((n + 3) * (n + 2) * (n + 1) * n);

      d2VdX2(0, 0) = cf(0) / 4 * (cf(7) * Vn2m2 - 3 * cf(8) * Vn2m);
      d2VdX2(0, 1) = cf(0) / 4 * (cf(7) * Wn2m2 - cf(8) * Wn2m);
      d2VdX2(1, 1) = cf(0) / 4 * (-cf(7) * Vn2m2 - cf(8) * Vn2m);

      d2WdX2(0, 0) = d2VdX2(0, 1);
      d2WdX2(0, 1) = d2VdX2(1, 1);
      d2WdX2(1, 1) = cf(0) / 4 * (-cf(7) * Wn2m2 - 3 * cf(8) * Wn2m);
    } else {
      if (m_in == 2) {
        kdel(1) = 1;
      } else {
        kdel(1) = 2;
      }

      cf(9) = sqrt((n + m + 4) * (n + m + 3) * (n + m + 2) * (n + m + 1));
      cf(10) = sqrt(2 / kdel(1) * (n - m + 4) * (n - m + 3) * (n - m + 2) *
                    (n - m + 1));
      d2VdX2(0, 0) =
          cf(0) / 4 * (cf(9) * Vn2m2 - 2 * cf(1) * Vn2m + cf(10) * Vn2m_2);
      d2VdX2(0, 1) = cf(0) / 4 * (cf(9) * Wn2m2 - cf(10) * Wn2m_2);
      d2VdX2(1, 1) =
          cf(0) / 4 * (-cf(9) * Vn2m2 - 2 * cf(1) * Vn2m - cf(10) * Vn2m_2);

      d2WdX2(0, 0) =
          cf(0) / 4 * (cf(9) * Wn2m2 - 2 * cf(1) * Wn2m + cf(10) * Wn2m_2);
      d2WdX2(0, 1) = cf(0) / 4 * (-cf(9) * Vn2m2 + cf(10) * Vn2m_2);
      d2WdX2(1, 1) =
          cf(0) / 4 * (-cf(9) * Wn2m2 - 2 * cf(1) * Wn2m - cf(10) * Wn2m_2);
    }
  }
  // === MAKE FULL MATRIX
  d2VdX2(1, 0) = d2VdX2(0, 1);
  d2VdX2(2, 0) = d2VdX2(0, 2);
  d2VdX2(2, 1) = d2VdX2(1, 2);
  d2WdX2(1, 0) = d2WdX2(0, 1);
  d2WdX2(2, 0) = d2WdX2(0, 2);
  d2WdX2(2, 1) = d2WdX2(1, 2);

  return std::tuple<Matrix3d, Matrix3d>(d2VdX2, d2WdX2);
}

Vector3d spharm_acc_ecr(int nmax, int mmax, const Vector3real &x_R_in,
                        double Re, double GMe, const MatrixXd &Cnm,
                        const MatrixXd &Snm) {
  Vector3d x_R = x_R_in.cast<double>();
  Vector3d accECR(0, 0, 0);
  MatrixXd V = MatrixXd::Zero(nmax + 3, mmax + 3);
  MatrixXd W = MatrixXd::Zero(nmax + 3, mmax + 3);
  V(0, 0) = Re / x_R.norm();

  Vector3d a, b;

  for (int n = 1; n < nmax + 3; ++n) {
    for (int m = 0; m < std::min(n + 1, mmax + 3); ++m) {
      if (n == (m + 1))
        std::tie(V(n, m), W(n, m)) = spharm_vwm1m(m, V(m, m), W(m, m), x_R, Re);
      else if (n == m)
        std::tie(V(n, m), W(n, m)) =
            spharm_vwmm(m, V(m - 1, m - 1), W(m - 1, m - 1), x_R, Re);
      else
        std::tie(V(n, m), W(n, m)) = spharm_vwnm(
            n, m, V(n - 1, m), V(n - 2, m), W(n - 1, m), W(n - 2, m), x_R, Re);

      if ((n > 2) && (m > 3)) {
        std::tie(a, b) = spharm_dvwdx(
            n - 2, m - 2, V(n - 1, m - 2), V(n - 1, m - 1), V(n - 1, m - 3),
            W(n - 1, m - 2), W(n - 1, m - 1), W(n - 1, m - 3), Re);
        accECR += (Cnm(n - 2, m - 2) * a + Snm(n - 2, m - 2) * b);
      } else if ((n > 2) && (m == 3)) {
        std::tie(a, b) = spharm_dvwdx(
            n - 2, m - 2, V(n - 1, m - 2), V(n - 1, m - 1), V(n - 1, m - 3),
            W(n - 1, m - 2), W(n - 1, m - 1), W(n - 1, m - 3), Re);
        accECR += (Cnm(n - 2, m - 2) * a + Snm(n - 2, m - 2) * b);
      } else if ((n >= 2) && (m == 2)) {
        std::tie(a, b) =
            spharm_dvwdx(n - 2, m - 2, V(n - 1, m - 2), V(n - 1, m - 1), 0,
                         W(n - 1, m - 2), W(n - 1, m - 1), 0, Re);
        accECR += (Cnm(n - 2, m - 2) * a + Snm(n - 2, m - 2) * b);
      }
    }
  }

  accECR *= (GMe / Re);
  return accECR;
}

}  // namespace lupnt