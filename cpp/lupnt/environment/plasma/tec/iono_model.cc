/**
 * @file tec/iono_model.cc
 * @brief Global selector for the electron-density backend used by compute_ne().
 *
 * @copyright Copyright (c) 2025
 */
#include "lupnt/environment/plasma/tec/iono_model.h"

namespace pecsim {

  // Default to GCPM so existing behavior is unchanged unless explicitly switched.
  static IonoModel g_iono_model = IonoModel::GCPM;

  void set_iono_model(IonoModel model) { g_iono_model = model; }

  IonoModel get_iono_model() { return g_iono_model; }

  std::string get_iono_model_str() {
    switch (g_iono_model) {
      case IonoModel::GCPM: return "GCPM";
      case IonoModel::NEQUICK_G: return "NeQuickG";
      case IonoModel::NEDM2020: return "NEDM2020";
    }
    return "GCPM";
  }

}  // namespace pecsim
