/**
 * @file tec/iono_model.h
 * @brief Global selector for the electron-density backend used by compute_ne().
 *
 * The ray tracer pulls electron density through a single choke point
 * (compute_ne). This selector chooses which backend that choke point
 * dispatches to. GCPM is the default so existing behavior is unchanged.
 *
 * NEQUICK_G is only functional when LuPNT is built with
 * -DLUPNT_ENABLE_NEQUICK=ON; otherwise selecting it and tracing raises a
 * runtime error from the stub (see nequick_stub.cc). The enum value always
 * exists so the Python/C++ ABI is stable regardless of build config.
 *
 * @copyright Copyright (c) 2025
 */
#pragma once

#include <array>
#include <string>

namespace pecsim {

  /** Electron-density backend used by compute_ne(). */
  enum class IonoModel {
    GCPM,       ///< Global Core Plasma Model (default; plasmasphere-capable).
    NEQUICK_G,  ///< Galileo NeQuick-G (requires -DLUPNT_ENABLE_NEQUICK=ON).
    NEDM2020    ///< Neustrelitz Electron Density Model 2020 (MIT; E+F, plasmasphere pending).
  };

  /** How the NeQuick Effective Ionisation Level (Az, sfu) is determined.
   *  These config types live in the core (not the gated module) so
   *  RayTraceConfig can carry them in both MIT and NeQuick builds. */
  enum class NeQuickAzMode {
    FROM_F107,  ///< Az = F10.7 from the IRI/GCPM pipeline (comparable to GCPM).
    EXPLICIT    ///< Az from az_sfu, or from ai[] evaluated at the point's MODIP.
  };

  /** Solar-activity driver for the NeQuick-G backend. */
  struct NeQuickSolarConfig {
    NeQuickAzMode mode = NeQuickAzMode::FROM_F107;
    double az_sfu = -1.0;  ///< Constant Az [sfu] when EXPLICIT and >= 0.
    std::array<double, 3> ai
        = {-1.0, 0.0, 0.0};  ///< ai0/ai1/ai2 (Az vs MODIP), used if az_sfu < 0.
  };

  /** Select the electron-density backend. Default is IonoModel::GCPM. */
  void set_iono_model(IonoModel model);

  /** @return the currently selected electron-density backend. */
  IonoModel get_iono_model();

  /** @return the name of the currently selected backend ("GCPM" | "NeQuickG"). */
  std::string get_iono_model_str();

}  // namespace pecsim
