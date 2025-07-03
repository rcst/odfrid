// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "adjust_with_caps.h"
#include "round_with_preserved_sum.h"
#include "uniform_simplex.h"
#include "math.h"
#include "likelihood.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]

// define global model parameters
arma::imat x;
arma::vec t;
arma::imat y;
arma::mat psi;
arma::mat phi;
double rho;

//' load - No. passengers on the bus immediatly after each stops
//'
//' @param x column vector of boardings and alightings
//' @return vector of passengers loadings immediatly after each stops
// [[Rcpp::export]]
arma::imat load(arma::imat& x) {
  arma::uword S = x.n_rows / 2.0;
  arma::imat w(S, x.n_cols);

  for(arma::uword c = 0; c < x.n_cols; ++c) {
    for(arma::uword i = 0; i<S; ++i) {
      if(i==0) w(i, c) = x(i, c) - x(S+i, c);
      else w(i, c) = w(i-1, c) + x(i, c) - x(S+i, c);
    }
  }

  return w;
}

//' routing_matrix - construct routing matrix from no. stops (S)
//'
//' @param S length-one integer denoting the number of stops
//' @return An (armadillo) integer matrix
// [[Rcpp::export]]
arma::imat routing_matrix(int s) {

  int m = ((s * s) - s) / 2;
  arma::imat A(2*s, m);

  // column index variable of A
  int j = 0;

  for(int i = 0; i < s; ++i) {
    for(int k = i + 1; k < s; ++k) {
      j = s * i - (i + 1) * (i + 2) / 2 + k;
      // std::cout << "i: " << i << ", j: " << j << ", k: " << k << std::endl;
      A(i, j) = 1.0;
    }
  }

  for(int i = s; i < 2*s; ++i) {
    for(int k = 0; k < (i + 1) - s - 1; ++k) {
      j = s * (k - 1) - (k + 1) * (k + 2) / 2 + i;
      // std::cout << "i: " << i << ", j: " << j << ", k: " << k << std::endl;
      A(i, j) = 1.0;
    }
  }

  return A;
}

//' ztoy - capped uniform simplex sampling 
//'
//' Helper function that generates an integer vector with constant total sum and varying caps
//'
//' @param z numeric vector of no. passengers approaching a single specific stop
//' @param v double no. alighters at that specific stop
//' @return An integer vector of same length as z whos values are all smaller-or-equal to z and it's sum is equal to v
// [[Rcpp::export]]
arma::ivec ztoy(arma::ivec z, double v) {
  std::vector<double> y = uniformSimplexSample(z.n_elem, v);
  std::vector<int> iy = roundWithPreservedSum(y);
  // iy = adjustWithCaps(iy, as<std::vector<int>>(z));
  iy = adjustWithCaps(iy, arma::conv_to<std::vector<int>>::from(z));
  arma::ivec r(iy.data(), iy.size(), true);
  return r;
}

arma::umat odform2sub(arma::umat& K, int j, int c) {
  arma::uvec C(j);
  C.fill(c);

  return join_rows(K(arma::span(0, j-1), j), C).st();
}

//' rod - Conditional Sampling of OD vectors
//' 
//' @param x a integer matrix whoes columns each contain boarding and alighting counts of 1 bus
//' journey
//' @return A named list of containing (1) the sampled OD vector (named y), (2) a corresponging vector (named z)
//' the log probability density from Markov chain transition probabilities (named lq)
// [[Rcpp::export]]
List rod(arma::imat& x) {
  arma::uword N = x.n_cols;
  arma::uword S = x.n_rows / 2;
  arma::uword D = S * (S-1) / 2;

  arma::imat u = x.rows(0, S-1);
  arma::imat v = x.rows(S, x.n_rows - 1);
  arma::imat y(D, N);
  arma::imat z(D, N);

  // tracking indices of vectorized matrices 
  arma::umat K(S, S);
  arma::uvec eids = arma::trimatl_ind(size(K), -1);
  K.elem(eids) = arma::linspace<arma::uvec>(0, D-1, D);
  K = K.st();

  arma::uword k =  0;
  arma::uword k_ijm1 = 0;

  for(arma::uword c = 0; c < N; ++c) {
    for(arma::uword j = 1; j < S; ++j) {
      // +++ SETTING Z +++ 
      for(arma::uword i = 0; i < j; ++i) {
        // k = ij_to_id(i, j, S);
        k = K(i, j);
        // k_ijm1= ij_to_id(i, j-1, S);
        k_ijm1 = K(i, j-1);

        if(i==j-1) z(k, c) = u(i, c);
        else z(k, c) = z(k_ijm1, c) - y(k_ijm1, c);
      }

      // +++ SETTING Y +++
      arma::uvec klist = arma::sub2ind(arma::size(y), odform2sub(K, j, c));
      y.elem(klist) = ztoy(z.elem(klist), v(j, c));
    }
  }

  // Markov Chain Transition Probabilities
  arma::imat w = load(x);

  arma::mat pi_1 = -1.0 * log_choose_mat(w.rows(arma::span(0, w.n_rows - 2)), v.rows(arma::span(1, v.n_rows - 1)));
  arma::mat pi_2 = log_choose_mat(z, y);

  return List::create(
      Named("y") = y,
      Named("z") = z,
      Named("q") = arma::sum(pi_1) + arma::sum(pi_2));
}

//' @param k covariance matrix for column psi_d
//' @param psi mapping-factor matrix
//' @param d the index of column of matrix psi to be sampled
// void ess_psi(arma::mat& K, arma::uword& d) {
//   arma::uword N = psi.n_cols;
//   arma::vec nu = mvnrnd(arma::zeros(N), K);
// 
//   double gamma = arma::randu();
// 
//   // here we need all model parameters
//   // what is a elegant way to avoid
//   // passing all parameters all the time?
//   // --> perhaps best not to make it a function
//   double log_c = 
//   
//   do {
//   } while ();
// 
// }

// [[Rcpp::export]]
arma::mat departure_times_covariance_matrix(const arma::vec& t, double sigma, double l) {
  arma::mat T1 = repmat(t, 1, t.n_elem);
  arma::mat T2 = repmat(t.t(), t.n_elem, 1);
  arma::mat K = std::pow(sigma, 2) * exp(square(T1 - T2) / (-2.0 * std::pow(l, 2)));      

  return K;
}
