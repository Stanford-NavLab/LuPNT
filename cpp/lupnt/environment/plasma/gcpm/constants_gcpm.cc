/**
 * @file constants.cpp
 * @author Keidai Iiyama
 * @brief  This file contains the implementation of the constants used in the
 * GCPM model
 * @version 0.1
 * @date 2025-02-14
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {
  double PN(int i, int j) { return origPN[i + j * 72]; }

  double PS(int i, int j) { return origPS[i + j * 72]; }
}  // namespace pecsim
