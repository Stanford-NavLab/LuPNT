#include "lupnt/dynamics/dynamics.h"

#include "lupnt/core/progress_bar.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {
  Ptr<IState> OrbitDynamics::PropagateState(const Ptr<IState> &state, Real t0, Real tf,
                                            MatXd *stm) {
    int repres = (int)state->GetStateType();
    assert(repres >= (int)OrbitStateRepres::CARTESIAN && repres <= (int)OrbitStateRepres::RTN
           && "Invalid OrbitState representation");

    OrbitState orbit_state = *std::static_pointer_cast<OrbitState>(state).get();
    Ptr<IState> state_new;
    if (stm == nullptr) {
      OrbitState orbit_state_new = PropagateState(orbit_state, t0, tf);
      state_new = std::make_shared<OrbitState>(orbit_state_new);
    } else {
      assert(stm->rows() == 6 && stm->cols() == 6 && "Invalid STM size");
      Mat6d stm_6;
      OrbitState orbit_state_new = PropagateState(orbit_state, t0, tf, &stm_6);
      *stm = stm_6;
      state_new = std::make_shared<OrbitState>(orbit_state_new);
    }
    return state_new;
  }

  VecX OrbitDynamics::Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm) {
    assert(x0.size() == 6 && "Invalid state size");
    Vec6 x0_6 = x0;
    Vec6 xf;
    if (stm == nullptr) {
      xf = Propagate(x0_6, t0, tf);
    } else {
      assert(stm->rows() == 6 && stm->cols() == 6 && "Invalid STM size");
      Mat6d stm_6;
      xf = Propagate(x0_6, t0, tf, &stm_6);
      *stm = stm_6;
    }
    return xf;
  }

  MatX6 OrbitDynamics::Propagate(const MatX6 &x0, Real t0, Real tf) {
    MatX6 xf = MatX6::Zero(x0.rows(), x0.cols());
    for (int i = 0; i < x0.cols(); i++) {
      Vec6 x0_i = x0.col(i);
      Vec6 xf_i = Propagate(x0_i, t0, tf);
      xf.row(i) = xf_i;
    }
    return xf;
  }

}  // namespace lupnt