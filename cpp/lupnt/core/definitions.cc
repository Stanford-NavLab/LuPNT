#include "lupnt/core/definitions.h"

#include "lupnt/core/error.h"

namespace lupnt {
  Real lupnt_epoch = 0.0;

  Real GetLupntEpoch() { return lupnt_epoch; }
  void SetLupntEpoch(Real epoch) { lupnt_epoch = epoch; }
}  // namespace lupnt
