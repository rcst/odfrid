#include "RcppArmadillo.h"
#include "globals.h"
#include "math.h"

using namespace arma;

//' L - log likelihood evaluation function
//'
//' Evalutes likelihood function for each OD vector given all current parameters.
//' Works for either multiple-trip OD vector matrix or single-column OD-vector-matrix.
//'
//' @param 
//' @return An vec object with as many entries as columns of y
vec log_likelihood(const imat& y, const imat& x, const mat& phi, const mat& psi) {
  int S = x.size() / 2;
  imat u = x(span(0, S-1), span::all);

  // columns -> trips
  // rows -> stops (S-2)
  mat G = phi * psi.t();
  mat lbd = matrix_softmax(G, rho);

  vec lq(y.n_cols);
  for(uword c = 0; c < y.n_cols; ++c)
    lq(c) = accu(lgamma(u + 1.0)) + accu((y.col(c) * log(lbd.col(c))) - lgamma(y.col(c) + 1.0));
  
  return lq;
}

vec log_likelihood(const imat& y, const imat& x) {
  return log_likelihood(y, x, phi, psi);
}

vec log_likelihood(const mat& psi) {
  return log_likelihood(y, x, phi, psi);
}
