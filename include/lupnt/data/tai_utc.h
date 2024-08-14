#include <lupnt/core/constants.h>

namespace lupnt {
  struct TaiUtcFileData {
    VecXd jd;
    VecXd tai_utc;
    VecXd mjd0;
    VecXd scale;
  };

  double GetTaiUtcDifference(double mjd);
}  // namespace lupnt
