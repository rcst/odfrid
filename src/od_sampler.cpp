#include "RcppArmadillo.h"
#include "odfrid/globals.h"
#include "odfrid/math.h"
#include "odfrid/uniform_simplex.h"
#include "odfrid/likelihood.h"

using namespace arma;

//' load - No. passengers on the bus immediatly after each stops
//'
//' @param x column vector of boardings and alightings
//' @return vector of passengers loadings immediatly after each stops
// [[Rcpp::export]]
arma::imat load(arma::imat& x) {
  uword S = x.n_rows / 2.0;
  imat w(S, x.n_cols);

  for(uword c = 0; c < x.n_cols; ++c) {
    for(uword i = 0; i<S; ++i) {
      if(i==0) w(i, c) = x(i, c) - x(S+i, c);
      else w(i, c) = w(i-1, c) + x(i, c) - x(S+i, c);
    }
  }

  return w;
}

//' ztoy - capped uniform simplex sampling 
//'
//' Helper function that generates an integer vector with constant total sum and varying caps
//'
//' @param z numeric vector of no. passengers approaching a single specific stop
//' @param v double no. alighters at that specific stop
//' @return An integer vector of same length as z whos values are all smaller-or-equal to z and it's sum is equal to v
ivec ztoy(ivec z, double v) {
  std::vector<double> y = uniformSimplexSample(z.n_elem, v);
  std::vector<int> iy = roundWithPreservedSum(y);
  // iy = adjustWithCaps(iy, as<std::vector<int>>(z));
  iy = adjustWithCaps(iy, conv_to<std::vector<int>>::from(z));
  ivec r(iy.data(), iy.size(), true);
  return r;
}

umat odform2sub(umat& K, int j, int c) {
  uvec C(j);
  C.fill(c);

  return join_rows(K(span(0, j-1), j), C).st();
}

//' rod - Conditional Sampling of OD vectors
//' 
//' @param x a integer matrix whoes columns each contain boarding and alighting counts of 1 bus
//' journey
//' @return A named list of containing (1) the sampled OD vector (named y), (2) a corresponging vector (named z)
//' the log probability density from Markov chain transition probabilities (named lq)
Rcpp::List rod(imat& x) {
  uword N = x.n_cols;
  uword S = x.n_rows / 2;
  uword D = S * (S-1) / 2;

  imat u = x.rows(0, S-1);
  imat v = x.rows(S, x.n_rows - 1);
  imat y(D, N);
  imat z(D, N);

  // tracking indices of vectorized matrices 
  umat K(S, S);
  uvec eids = trimatl_ind(size(K), -1);
  K.elem(eids) = linspace<uvec>(0, D-1, D);
  K = K.st();

  uword k =  0;
  uword k_ijm1 = 0;

  for(uword c = 0; c < N; ++c) {
    for(uword j = 1; j < S; ++j) {
      // +++ SETTING Z +++ 
      for(uword i = 0; i < j; ++i) {
        // k = ij_to_id(i, j, S);
        k = K(i, j);
        // k_ijm1= ij_to_id(i, j-1, S);
        k_ijm1 = K(i, j-1);

        if(i==j-1) z(k, c) = u(i, c);
        else z(k, c) = z(k_ijm1, c) - y(k_ijm1, c);
      }

      // +++ SETTING Y +++
      uvec klist = sub2ind(size(y), odform2sub(K, j, c));
      y.elem(klist) = ztoy(z.elem(klist), v(j, c));
    }
  }

  // Markov Chain Transition Probabilities
  imat w = load(x);

  mat pi_1 = -1.0 * log_choose_mat(w.rows(span(0, w.n_rows - 2)), v.rows(span(1, v.n_rows - 1)));
  mat pi_2 = log_choose_mat(z, y);

  return Rcpp::List::create(
      Rcpp::Named("y") = y,
      Rcpp::Named("z") = z,
      Rcpp::Named("q") = sum(pi_1) + sum(pi_2));
}

