/**
 * @file tec/nequick_stub.cc
 * @brief Fallback definition of compute_ne_nequick when the NeQuick-G backend
 *        is not compiled in.
 *
 * This translation unit is ALWAYS compiled. When LuPNT is built with
 * -DLUPNT_ENABLE_NEQUICK=ON, LUPNT_HAS_NEQUICK is defined and the real
 * definition in plasma/nequick/ is used instead (this file then compiles to
 * nothing). When the NeQuick module is absent (e.g. an MIT release build),
 * this stub provides a single definition that raises a clear error if the
 * NeQuick backend is selected.
 *
 * Keeping the fallback here — rather than behind #ifdefs scattered through
 * raytrace.cc — lets the dispatcher in compute_ne() stay backend-agnostic.
 *
 * @copyright Copyright (c) 2025
 */
#ifndef LUPNT_HAS_NEQUICK

#  include <stdexcept>

#  include "lupnt/environment/plasma/tec/raytrace.h"

namespace pecsim {

  double compute_ne_nequick(double /*t_j2000*/, const Vec3d& /*pos_geo*/, RayTraceConfig /*config*/,
                            bool /*debug*/) {
    throw std::runtime_error(
        "NeQuick-G electron-density backend is not compiled in. "
        "Reconfigure LuPNT with -DLUPNT_ENABLE_NEQUICK=ON (requires the "
        "plasma/nequick/ module) to use IonoModel::NEQUICK_G.");
  }

}  // namespace pecsim

#endif  // !LUPNT_HAS_NEQUICK
