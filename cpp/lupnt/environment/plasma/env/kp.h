/**
 * @file env/kp.h
 * @author Keidai Iiyama
 * @brief This file contains the interface for retrieving the Kp index.
 * @version 0.1
 * @date 2025-02-17
 */

#pragma once

#include "lupnt/environment/plasma/env/time_utils.h"

namespace pecsim {

  /**
   * @brief Get the Kp index for a given datetime.
   * The Kp index is a measure of geomagnetic activity.
   *
   * @param datetime The date and time for which to retrieve the Kp index.
   * @return double The Kp index value.
   */
  double get_kp_index(DateTime datetime);

}  // namespace pecsim
