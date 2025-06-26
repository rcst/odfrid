// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "adjust_with_caps.h"
#include "round_with_preserved_sum.h"
#include "uniform_simplex.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// //' ij_to_id - OD-vector(s) index mapping function
// //' 
// //' Help function to flatten indexing of OD-vectors and other vectors of the
// //' same structure
// //'
// //' @param i (zero-based) index of starting-stop
// //' @param j (zero-based) index of ending-stop
// //' @param S total no. stops
// //' @return A "flattened" one-dimensional index that enumerates the vector
// //' entries
// //' (3) a vector of corresponding Markov chain transition probabilities
// arma::uword ij_to_id(arma::uword i, arma::uword j, arma::uword S) {
//   return (2*i*S - i*i + 2*j - 3*i - 2) / 2;
// }

arma::vec log_choose_vec(const arma::ivec& n, const arma::ivec& k) {
  if (n.n_elem != k.n_elem) {
    stop("Vectors 'n' and 'k' must have the same length.");
  }

  arma::vec n_d = arma::conv_to<arma::vec>::from(n);
  arma::vec k_d = arma::conv_to<arma::vec>::from(k);
  arma::vec nk_d = n_d - k_d;

  arma::vec result(n_d.n_elem, arma::fill::none);

  // Logical masks
  arma::uvec is_valid = (k >= 0) && (k <= n);
  arma::uvec is_equal = (n == k);
  arma::uvec is_zero  = (k == 0);
  arma::uvec compute  = is_valid && is_equal == 0 && is_zero == 0;

  // Set results
  result.fill(-arma::datum::inf);  // default for invalid cases

  // n == k → log(choose(n, n)) = 0
  result.elem(find(is_equal)).zeros();

  // k == 0 → log(choose(n, 0)) = 0
  result.elem(find(is_zero)).zeros();

  // General case
  arma::uvec idx = find(compute);
  if (!idx.is_empty()) {
    arma::vec n_sub  = n_d.elem(idx);
    arma::vec k_sub  = k_d.elem(idx);
    arma::vec nk_sub = nk_d.elem(idx);

    arma::vec logc = arma::lgamma(n_sub + 1.0) - arma::lgamma(k_sub + 1.0) - arma::lgamma(nk_sub + 1.0);
    result.elem(idx) = logc;
  }

  return result;
}

// [[Rcpp::export]]
arma::mat log_choose_mat(const arma::imat& N, const arma::imat& K) {
  if (size(N) != size(K)) {
    stop("N and K must be the same size.");
  }

  int n_rows = N.n_rows;
  int n_cols = N.n_cols;
  arma::mat result(n_rows, n_cols, arma::fill::none);

  for (int col = 0; col < n_cols; ++col) {
    arma::ivec n_col = N.col(col);
    arma::ivec k_col = K.col(col);

    // Use the optimized log_choose_vector function
    arma::vec logc = log_choose_vec(n_col, k_col);

    result.col(col) = logc;
  }

  return result;
}

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

  NumericVector pi(S);

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

// // [[Rcpp::export]]
// NumericVector softmax(const NumericVector &x) {
//   NumericVector ex;
//   ex = exp(x);
//   return ex / (1.0 * sum(ex));
// }

// [[Rcpp::export]]
arma::vec softmax(const arma::vec &x) {
  arma::vec y(x.size());
  y = exp(x);
  y(y.size() - 1) = 1; 
  y = y / (1.0 * sum(y));
  return y;
}


//' Compute softmax column-wise with fixed numerator = 1 for last entry of each segment
//' 
//' NOTE: Lambda should be returned so that the entry from second to last to last is always 1
//' 
//' @param G: M x N matrix (G = Phi %*% t(Psi))
//' @param rho: temperature parameter
//' @return: M x N matrix Lambda of softmax values
arma::mat matrix_softmax(const arma::mat& G, double rho) {
  // G has one row less than lambda - ie, G_S-1
  int M = G.n_rows;
  int N = G.n_cols;
  int S = (1.0 + sqrt(1.0 + (8.0 * (M+1.0)))) / 2.0;

  arma::ivec segment_lengths(S-1);
  segment_lengths = rev(seq(2, S-1));

  arma::mat Lambda(M+1, N, arma::fill::ones);
  int start = 0;

  for (int i = 0; i < S-2; ++i) {
    int len = segment_lengths(i);
    int end = start + len;

    for (int col = 0; col < N; ++col) {
      arma::vec g_block = G.submat(start, col, end - 1, col);
      arma::vec scaled = rho * g_block;

      arma::vec exp_vals(len);
      exp_vals(arma::span(0, len - 2)) = exp(scaled(arma::span(0, len - 2)) - max(scaled));  // stabilization
      exp_vals(len - 1) = 1.0;

      double denom = 1.0 + accu(exp_vals);
      for (int k = 0; k < len; ++k) {
        Lambda(start + k, col) = exp_vals(k) / denom;
      }
    }

    start += len;
  }

  return Lambda;
}

// [[Rcpp::export]]
double log_likelihood(NumericVector y, NumericVector x, NumericMatrix phi, NumericMatrix psi, double rho) {

  // perhaps best is to have this function
  // evaluate the likelihood of a single y^n
  // instead of all trips together -> needed for metropolis-hastings
  //
  // then have another function which computes the overall likelihood
  // as the product of all p(y^n | lambda^n)

  int S = x.size() / 2;
  // alternative arma::uvec (for unsigned int vectors)
  arma::vec ay = as<arma::vec>(y);
  // arma::vec ax = as<arma::vec>(x);
  arma::vec u = as<arma::vec>(x)(arma::span(0, S-1));

  arma::mat aphi = as<arma::mat>(phi);
  arma::mat apsi = as<arma::mat>(psi);

  // columns -> trips
  // rows -> stops (S-2)
  arma::mat G = aphi * apsi.t();
  arma::mat lambda = matrix_softmax(G, rho);

  double logl = accu(lgamma(u + 1.0));
  // logl = accu(lgamma(x(arma::span(0, S-1)))) + accu(y * log(lambda) - lgamma(y + 1.0));
  // logl = accu(lgamma(u + 1.0)) + accu(y * log(lambda) - lgamma(y + 1.0));
  for(arma::uword c = 0; c < lambda.n_cols; ++c) {
    logl += accu(ay.col(c) * log(lambda.col(c)) - lgamma(ay.col(c) + 1.0));
  }
  
  // Eq (6.11)
  // for each bus trip
  // double log_sum_u = sum(lfactorial(x[seq(0, S-1)]));
  // double log_sum_u = accu(lgamma(x(arma::span(0, S-1))));
  // double log_sum_mult = sum(log(pow(lambda, y)) - lfactorial(y));
  // double log_sum_mult = accu(y * log(lambda) - lgamma(y + 1.0)); // missing power of y_ij!

  // arma::vec likelihood(lambda.n_cols);
  // for(arma::uword n = 0; n < lambda.n_cols; ++n) {
  //   log_sum_u = sum(lfactorial(x[seq(0, S-1)]));
  //   log_sum_mult = sum(log(pow(lambda[,n], y)) - lfactorial(y));
  // }


  // return log_sum_u + log_sum_mult;
  return logl;
}
