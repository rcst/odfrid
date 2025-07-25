#include "RcppArmadillo.h"
#include "odfrid/globals.h"
#include "odfrid/math.h"

using namespace arma;

//' L - log likelihood evaluation function
//'
//' Evalutes likelihood function for each OD vector given all current parameters.
//' Works for either multiple-trip OD vector matrix or single-column OD-vector-matrix.
//'
//' @param 
//' @return An vec object with as many entries as columns of y
double log_likelihood(const umat& y, const umat& x, const mat& phi, const mat& psi, const double& rho, const umat& A) {
  uword S = x.n_rows / 2;
  umat u = x.rows(0, S-2);

  // columns -> trips
  // rows -> stops (S-2)
  mat G = phi * psi.t();
  mat lbd = matrix_softmax(G, rho);

  umat x_check = A * y;

  vec lq(y.n_cols);
  lq.fill(-datum::inf);

  urowvec lq_1 = sum(lgamma(u + 1.0));
  rowvec lq_2 = sum((y % log(lbd)) - lgamma(y));

  // according to the size of y 
  // to enable single-trip-likelihood
  for(uword c = 0; c < y.n_cols; ++c)
    if(all(x_check.col(c) == x.col(c)))
      lq(c) = lq_1(c) + lq_2(c);

  // parallelize!
  return accu(lq);
}

double log_likelihood(const uvec& y, const uvec& x, const uword n) {
  return log_likelihood(y, x, phi, psi.row(n), rho, A);
}

double log_likelihood(const umat& y, const umat& x) {
  return log_likelihood(y, x, phi, psi, rho, A);
}

double log_likelihood(const mat& psi, const uword d) {
  return log_likelihood(y, x, phi, psi, rho, A);
}

double log_likelihood(const uword d, const mat& phi) {
  return log_likelihood(y, x, phi, psi, rho, A);
}

double log_likelihood(const double& rho) {
  return log_likelihood(y, x, phi, psi, rho, A);
}

double log_likelihood() {
  return log_likelihood(y, x, phi, psi, rho, A);
}
