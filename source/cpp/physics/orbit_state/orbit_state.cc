#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  std::ostream &operator<<(std::ostream &os, const OrbitStateRepres &repres) {
    switch (repres) {
      case OrbitStateRepres::CARTESIAN: os << "CARTESIAN"; break;
      case OrbitStateRepres::CLASSICAL_OE: os << "CLASSICAL_OE"; break;
      case OrbitStateRepres::QUASI_NONSINGULAR_OE: os << "QUASI_NONSINGULAR_OE"; break;
      case OrbitStateRepres::SINGULAR_ROE: os << "SINGULAR_ROE"; break;
      case OrbitStateRepres::NONSINGULAR_OE: os << "NONSINGULAR_OE"; break;
      case OrbitStateRepres::EQUINOCTIAL_OE: os << "EQUINOCTIAL_OE"; break;
      case OrbitStateRepres::DELAUNAY_OE: os << "DELAUNAY_OE"; break;
      case OrbitStateRepres::ABSOLUTE_RELATIVE_SEPARATOR:
        os << "ABSOLUTE_RELATIVE_SEPARATOR";
        break;
      case OrbitStateRepres::RTN: os << "RTN"; break;
      case OrbitStateRepres::QUASINONSINGULAR_ROE: os << "QUASINONSINGULAR_ROE"; break;
    }
    return os;
  }

  // ****************************************************************************
  // OrbitState
  // ****************************************************************************

  OrbitState::OrbitState(const Vec6 &x, Frame coord, OrbitStateRepres repres,
                         const std::array<const char *, 6> &names,
                         const std::array<const char *, 6> &units)
      : x_(x), frame_(coord), repres_(repres), names_(names), units_(units) {}

  // Overrides
  int OrbitState::GetSize() const { return kOrbitStateSize; }
  VecX OrbitState::GetVec() const { return x_; }
  void OrbitState::SetVec(const VecX &x) { x_ = x; }
  Real OrbitState::GetValue(int i) const { return x_(i); }
  void OrbitState::SetValue(int idx, Real val) { x_(idx) = val; }
  StateType OrbitState::GetStateType() const { return static_cast<StateType>(repres_); }

  Vec6 OrbitState::GetVec6() const { return x_; }
  Frame OrbitState::GetFrame() const { return frame_; }
  std::array<const char *, 6> OrbitState::GetUnits() const { return units_; }
  std::array<const char *, 6> OrbitState::GetNames() const { return names_; }
  OrbitStateRepres OrbitState::GetOrbitStateRepres() const { return repres_; }

  void OrbitState::SetOrbitStateRepres(const OrbitStateRepres rep) { repres_ = rep; }
  void OrbitState::SetCoordSystem(Frame frame) { frame_ = frame; }

  Real OrbitState::operator()(int idx) const { return x_(idx); }
  std::ostream &OrbitState::operator<<(std::ostream &os) const {
    os << "<OrbitState(" << x_.transpose() << ", " << frame_ << ", " << repres_ << ")>";
    return os;
  }

  // ****************************************************************************
  // CartesianOrbitState
  // ****************************************************************************

}  // namespace lupnt
