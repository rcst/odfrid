#include "Rcpp.h"

using namespace Rcpp;

// Generate N values uniformly distributed on the simplex (sum = 1),
// then optionally scale and shift each by C and A respectively.
// [[Rcpp::export]]
std::vector<double> uniformSimplexSample(int N, double C = 1.0, double A = 0.0) {
  // std::random_device rd;
  // std::mt19937 gen(rd()); // Mersenne Twister RNG
  // std::uniform_real_distribution<> uniform(0.0, 1.0);

  std::vector<double> u(N);
  for (int i = 0; i < N; ++i) {
    double v = (double)runif(1, 0, 1)[0];
    u[i] = -std::log(v); // Exponential(1) sample
  }

  double sumU = std::accumulate(u.begin(), u.end(), 0.0);

  std::vector<double> p(N);
  for (int i = 0; i < N; ++i) {
    p[i] = C * (u[i] / sumU) + A;
  }

  return p;
}
