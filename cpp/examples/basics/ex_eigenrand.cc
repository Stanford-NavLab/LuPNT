// Draws samples from a correlated multivariate normal distribution with
// EigenRand and plots the sample cloud.
#include <lupnt/lupnt.h>
using namespace lupnt;

int main() {
  RandomEngine::SetSeed(1234);

  Vec2d mean{0.0, 1.0};

  double sig1 = 1.0;
  double sig2 = 3.0;
  double rho = 0.8;
  Mat2d cov{{sig1 * sig1, rho * sig1 * sig2}, {rho * sig1 * sig2, sig2 * sig2}};

  MatX2d samples = SampleMvNormal(mean, cov, 1000).cast<double>();

  matplot::plot(samples.col(0).eval(), samples.col(1).eval(), "o");
  matplot::show();

  return 0;
}
