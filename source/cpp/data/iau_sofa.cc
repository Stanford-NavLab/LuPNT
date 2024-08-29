#include "lupnt/data/iau_sofa.h"

#include <lupnt/core/constants.h>
#include <lupnt/core/file.h>
#include <lupnt/numerics/interpolation.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

namespace lupnt {

  Ptr<IauSofaFileData> iau_sofa;
  std::mutex iau_sofa_mutex;

  void LoadIauSofaFileData(const std::filesystem::path& filepath) {
    std::lock_guard<std::mutex> lock(iau_sofa_mutex);
    if (iau_sofa) return;  // Data already loaded

    size_t n_lines = CountLines(filepath.string());
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    iau_sofa = MakePtr<IauSofaFileData>();
    iau_sofa->jd_tt.resize(n_lines);
    iau_sofa->X.resize(n_lines);
    iau_sofa->Y.resize(n_lines);
    iau_sofa->s.resize(n_lines);

    size_t row = 0;
    std::string line;
    double jd_tt, X, Y, s;
    // Read data lines
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      iss >> jd_tt >> X >> Y >> s;
      iau_sofa->jd_tt(row) = jd_tt;
      iau_sofa->X(row) = X;
      iau_sofa->Y(row) = Y;
      iau_sofa->s(row) = s;
      ++row;
    }
    file.close();
    return;
  }

  IauSofaData GetIauSofaData(Real jd_tt) {
    if (!iau_sofa) LoadIauSofaFileData(GetFilePath(IAU_SOFA_FILENAME));

    IauSofaData data;

    size_t index;
    if (jd_tt < iau_sofa->jd_tt(0) || jd_tt > iau_sofa->jd_tt(iau_sofa->jd_tt.size() - 1)) {
      index = (jd_tt < iau_sofa->jd_tt(0)) ? 0 : iau_sofa->jd_tt.size() - 1;
      data.X = iau_sofa->X(index);
      data.Y = iau_sofa->Y(index);
      data.s = iau_sofa->s(index);
    }

    int order = 9;
    LagrangeInterpolator interp(iau_sofa->jd_tt, jd_tt.val(), order);
    data.X = interp.Interpolate(iau_sofa->X);
    data.Y = interp.Interpolate(iau_sofa->Y);
    data.s = interp.Interpolate(iau_sofa->s);
    return data;
  }

}  // namespace lupnt
