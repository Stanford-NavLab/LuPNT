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

  enum class StateType : int;

  class IState {
  public:
    virtual ~IState() = default;
    virtual int GetSize() const = 0;
    virtual VecX GetVec() const = 0;
    virtual void SetVec(const VecX &x) = 0;
    virtual Real GetValue(int idx) const = 0;
    virtual void SetValue(int idx, Real val) = 0;
    virtual StateType GetStateType() const = 0;
  };

}  // namespace lupnt
