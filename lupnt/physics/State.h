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

#include <vector>



namespace lupnt {

/**
 * @brief Interface for States
 *
 */
class IState {
 public:
  virtual ~IState() = default;
  virtual inline int GetStateSize() const = 0;
  virtual inline real GetValue(int idx) const = 0;
};

/**
 * @brief Stack of multiple state types (example: orbit and clock)
 *
 */
class JointState {
 private:
  std::vector<IState*> state_vec_;
  int state_vec_size_ = 0;

 public:
  JointState(){};
  JointState(std::vector<IState*> state_vec) {
    int state_vec_size = 0;
    for (int i = 0; state_vec.size(); i++) {
      state_vec_.push_back(state_vec[i]);
      state_vec_size += state_vec_[i]->GetStateSize();
    }
    state_vec_size_ = state_vec_size;
  };

  void PushBackState(IState* state) { state_vec_.push_back(state); };

  std::vector<IState*> GetJointState() { return state_vec_; };

  VectorXreal GetJointStateValue() {
    VectorXreal advec(state_vec_size_);
    int cur_idx = 0;
    for (int i = 0; state_vec_.size(); i++) {
      for (int j = 0; j < state_vec_[i]->GetStateSize(); j++) {
        advec(i) = state_vec_[i]->GetValue(j);
      }
    }
    return advec;
  };
};

}  // namespace lupnt
