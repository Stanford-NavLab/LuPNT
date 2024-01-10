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

namespace lupnt {

struct Body {
  std::string name;
  double mu;
  double R;
  NaifId id;
  bool sphericalHarmonics;
  bool normalized;
  int n_max;
  int m_max;
  MatrixXd Cnm;
  MatrixXd Snm;

  static Body Moon(int n_max = 0, int m_max = 0);
  static Body Earth(int n_max = 0, int m_max = 0);
};

struct BodyData {
  std::string name;
  double GM;
  double R;

  std::string filepath;
  int headerlines;
  std::string delimiter;
};

std::vector<std::string> split(const std::string &s, char delimiter);

void ReadData(const std::string &filepath, int N, int headerlines,
              const std::string &delimiter, std::vector<int> &idN,
              std::vector<int> &idM, std::vector<double> &C,
              std::vector<double> &S, std::vector<double> &sigC,
              std::vector<double> &sigS);
              
BodyData GetBodyData(const NaifId bodyID);

std::tuple<MatrixXd, MatrixXd> LoadGravityCoefficients(BodyData bd, int nmax);

}  // namespace lupnt