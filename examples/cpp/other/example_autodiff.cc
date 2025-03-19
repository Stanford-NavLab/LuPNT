#include <lupnt/lupnt.h>

#include <iostream>
using namespace lupnt;
using namespace std;

int main() {
  Vec3 a(1.42312, 2.213453254325, 3.5342543543231245231);
  Vec3 b(4, 5, 6);
  Vec3 c(7, 8, 9);
  Mat3 m;
  m << a, b, c;
  cout << "Clean  " << endl << m.format(FMT_CLEAN) << endl << endl;
  cout << "Compact" << endl << m.format(FMT_COMPACT) << endl << endl;
  cout << "Heavy  " << endl << m.format(FMT_HEAVY) << endl << endl;
  return 0;
}
