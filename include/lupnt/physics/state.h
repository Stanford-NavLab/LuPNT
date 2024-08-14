/**
 * @file state.h
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

namespace lupnt {

  /**
   * @brief Interface for States
   *
   */
  class IState {
  public:
    virtual ~IState() = default;
    virtual inline int GetSize() const = 0;
    virtual inline Real GetValue(int idx) const = 0;
    virtual inline void SetValue(Real val, int idx) = 0;
  };

}  // namespace lupnt
