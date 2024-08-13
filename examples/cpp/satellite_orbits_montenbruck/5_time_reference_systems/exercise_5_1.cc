#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;

int main() {
  Real t_utc = GregorianToTime(1999, 3, 4, 0, 0, 0);
  Real t_tai = ConvertTime(t_utc, TimeSys::UTC, TimeSys::TAI);
}
