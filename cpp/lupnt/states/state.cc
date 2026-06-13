#include "lupnt/states/state.h"

#include <sstream>

#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"

namespace lupnt {

  State::State() : VecX() {}

  State::State(int n) : VecX(n) { SetNamesAndUnits(n); }

  State::State(int rows, int cols) : VecX(rows) {
    LUPNT_CHECK(cols == 1, "State must have 1 column", "State");
    SetNamesAndUnits(rows);
  }

  State::State(const VecX& vec, const std::vector<std::string>& names,
               const std::vector<std::string>& units)
      : VecX(vec) {
    names_ = names;
    units_ = units;
  }

  State::State(const State& other) : VecX(other) {
    type_ = other.type_;
    names_ = other.names_;
    units_ = other.units_;
    frame_ = other.frame_;
  }

  State& State::operator=(const State& other) {
    VecX::operator=(other);

    type_ = other.type_;
    names_ = other.names_;
    units_ = other.units_;
    frame_ = other.frame_;
    return *this;
  }

  StateType State::GetType() const { return type_; }
  std::string State::GetName() const { return name_; }
  std::vector<std::string> State::GetNames() const { return names_; }
  std::vector<std::string> State::GetUnits() const { return units_; }
  Frame State::GetFrame() const { return frame_; }

  void State::SetType(StateType type) { type_ = type; }
  void State::SetName(const std::string& name) { name_ = name; }
  void State::SetNames(const std::vector<std::string>& names) { names_ = names; }
  void State::SetUnits(const std::vector<std::string>& units) { units_ = units; }
  void State::SetFrame(Frame frame) { frame_ = frame; }

  std::string State::ToString() const {
    std::ostringstream oss;
    oss << type_ << "(" << this->transpose().format(FMT_COMPACT) << ", " << enum_name(frame_)
        << ")";
    return oss.str();
  }

  void State::SetNamesAndUnits(int n) {
    names_.resize(n);
    units_.resize(n);
    for (int i = 0; i < n; ++i) {
      names_[i] = std::to_string(i);
      units_[i] = "none";
    }
  }

}  // namespace lupnt
