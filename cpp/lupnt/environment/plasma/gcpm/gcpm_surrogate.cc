/**
 * @file gcpm/gcpm_surrogate.cc
 * @brief Loader + 7-D interpolation for the GCPM surrogate (see gcpm_surrogate.h).
 *
 * @copyright Copyright (c) 2026
 */
#include "lupnt/environment/plasma/gcpm/gcpm_surrogate.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "lupnt/environment/plasma/core/user_filepath.h"

namespace pecsim {

  namespace {
    // Bracket value v on ascending axis ax; returns (i0,i1,w) for linear interp.
    // Periodic axes wrap with period P (values assumed within one period).
    void bracket(const std::vector<double>& ax, bool per, double P, double v, int& i0, int& i1,
                 double& w) {
      int n = static_cast<int>(ax.size());
      if (n == 1) {
        i0 = i1 = 0;
        w = 0.0;
        return;
      }
      if (per) {
        v = std::fmod(v - ax[0], P);
        if (v < 0) v += P;
        v += ax[0];
        if (v >= ax[n - 1]) {
          i0 = n - 1;
          i1 = 0;
          w = (v - ax[n - 1]) / (ax[0] + P - ax[n - 1]);
          return;
        }
      } else {
        if (v <= ax[0]) {
          i0 = 0;
          i1 = 0;
          w = 0.0;
          return;
        }
        if (v >= ax[n - 1]) {
          i0 = n - 1;
          i1 = n - 1;
          w = 0.0;
          return;
        }
      }
      int i = 0;
      while (i < n - 1 && ax[i + 1] <= v) ++i;
      i0 = i;
      i1 = i + 1;
      w = (v - ax[i]) / (ax[i + 1] - ax[i]);
    }
  }  // namespace

  bool GcpmSurrogate::load(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    auto rd_i32 = [&] {
      int32_t x = 0;
      f.read(reinterpret_cast<char*>(&x), 4);
      return x;
    };
    auto rd_f64 = [&] {
      double x = 0;
      f.read(reinterpret_cast<char*>(&x), 8);
      return x;
    };

    if (rd_i32() != 0x47435053) return false;  // magic "GCPS"
    int ndim = rd_i32();
    if (ndim != 7) return false;
    axes_.resize(ndim);
    periodic_.assign(ndim, false);
    period_.assign(ndim, 0.0);
    long ntot = 1;
    for (int a = 0; a < ndim; ++a) {
      int len = rd_i32();
      periodic_[a] = rd_i32() != 0;
      period_[a] = rd_f64();
      axes_[a].resize(len);
      for (int k = 0; k < len; ++k) axes_[a][k] = rd_f64();
      ntot *= len;
    }
    // C-order strides.
    stride_.assign(ndim, 1);
    for (int a = ndim - 2; a >= 0; --a) stride_[a] = stride_[a + 1] * axes_[a + 1].size();
    table_.resize(ntot);
    f.read(reinterpret_cast<char*>(table_.data()), ntot * sizeof(float));
    loaded_ = static_cast<bool>(f);
    return loaded_;
  }

  double GcpmSurrogate::eval(double r, double lam_m, double mlt, double ut, double kp, double r12,
                             double doy) const {
    double q[7] = {r, lam_m, mlt, ut, kp, r12, doy};
    int i0[7], i1[7];
    double w[7];
    for (int a = 0; a < 7; ++a)
      bracket(axes_[a], periodic_[a], period_[a], q[a], i0[a], i1[a], w[a]);
    double acc = 0.0;
    for (int c = 0; c < 128; ++c) {
      double wt = 1.0;
      long idx = 0;
      for (int a = 0; a < 7; ++a) {
        int b = (c >> a) & 1;
        wt *= b ? w[a] : (1.0 - w[a]);
        idx += (b ? i1[a] : i0[a]) * stride_[a];
      }
      if (wt != 0.0) acc += wt * table_[idx];
    }
    return std::exp(acc);
  }

  const GcpmSurrogate& GcpmSurrogate::instance() {
    static GcpmSurrogate g = [] {
      GcpmSurrogate s;
      const char* env = std::getenv("LUPNT_GCPM_SURROGATE");
      std::string path = env ? std::string(env) : get_base_path() + "/gcpm_surrogate.bin";
      if (!s.load(path)) {
        std::fprintf(stderr, "[gcpm_surrogate] table not found at %s; falling back to full GCPM.\n",
                     path.c_str());
      }
      return s;
    }();
    return g;
  }

}  // namespace pecsim
