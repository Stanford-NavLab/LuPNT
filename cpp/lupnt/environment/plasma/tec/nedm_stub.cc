/**
 * @file tec/nedm_stub.cc
 * @brief Fallback definition of compute_ne_nedm when the NEDM2020 backend
 *        is not compiled in.
 *
 * This translation unit is ALWAYS compiled. When LuPNT is built with
 * -DLUPNT_ENABLE_NEDM=ON, LUPNT_HAS_NEDM is defined and the real definition
 * in plasma/nedm/ is used instead (this file then compiles to nothing). When
 * the NEDM module is absent (e.g. an MIT release build), this stub provides a
 * single definition that raises a clear error if the NEDM backend is selected.
 *
 * Keeping the fallback here — rather than behind #ifdefs scattered through
 * raytrace.cc — lets the dispatcher in compute_ne() stay backend-agnostic.
 *
 * @copyright Copyright (c) 2025
 */
#ifndef LUPNT_HAS_NEDM

#  include <stdexcept>

#  include "lupnt/environment/plasma/tec/raytrace.h"

namespace pecsim {

  double compute_ne_nedm(double /*t_j2000*/, const Vec3d& /*pos_geo*/, RayTraceConfig /*config*/,
                         bool /*debug*/) {
    throw std::runtime_error(
        "NEDM2020 electron-density backend is not compiled in. "
        "Reconfigure LuPNT with -DLUPNT_ENABLE_NEDM=ON (requires the "
        "plasma/nedm/ module) to use IonoModel::NEDM2020.");
  }

}  // namespace pecsim

#endif  // !LUPNT_HAS_NEDM
