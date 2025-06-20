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

//' ij_to_id - OD-vector(s) index mapping function
//' 
//' Help function to flatten indexing of OD-vectors and other vectors of the
//' same structure
//'
//' @param i (zero-based) index of starting-stop
//' @param j (zero-based) index of ending-stop
//' @param S total no. stops
//' @return A "flattened" one-dimensional index that enumerates the vector
//' entries
//' (3) a vector of corresponding Markov chain transition probabilities
int ij_to_id(int i, int j, int S) {
  return (2*i*S - i*i + 2*j - 3*i - 2) / 2;
}

//' i_to_id - subsetting helper function 
//'
//' Helper function to find all indices from a vector with boarding-alighting-index-structure
//'
//' @param i integer denoting the ith stop
//' @param S integer denoting the total no. stops
//' @return An (armadilllo) unsigned integeger vector of parameter indexes that correspond to boardings at the ith stop
// [[Rcpp::export]]
arma::umat i_to_id(unsigned int i, unsigned int N, unsigned int S){
  // first col - i
  // second col - n
  arma::umat ids(S-i-1, 2);
  for(unsigned int k = 0; k < (S-i-1); ++k) {
    ids(k, 0) = ij_to_id(i, i+k+1, S);
    ids(k, 1) = N;
  }

  return ids.t();
}

// arma::uvec i_to_id(unsigned int i, unsigned int S) {
//   arma::uvec ids(S-i-1);
//   for(unsigned int k = 0; k < (S-i-1); ++k)
//     ids(k) = ij_to_id(i, i+k+1, S);
// 
//   return ids;
// }

//' load - No. passengers on the bus immediatly after each stops
//'
//' @param x column vector of boardings and alightings
//' @return vector of passengers loadings immediatly after each stops
// [[Rcpp::export]]
IntegerVector load(IntegerVector x) {
  int S = x.size() / 2.0;
  // std::cout << "S: " << S << std::endl;
  IntegerVector w(S);

  for(int i = 0; i<S; ++i) {
    // std::cout << x(i) << " " << x(S+i) << std::endl;
    if(i==0) w(i) = x(i) - x(S+i);
    else w(i) = w(i-1) + x(i) - x(S+i);
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
IntegerVector ztoy(IntegerVector z, double v) {
  std::vector<double> y = uniformSimplexSample(z.length(), v);
  std::vector<int> iy = roundWithPreservedSum(y);
  iy = adjustWithCaps(iy, as<std::vector<int>>(z));
  return wrap(iy);
}

//' rod - Conditional Sampling of OD vectors
//' 
//' @param x a column vector containing boarding and alighting counts of 1 bus
//' journey
//' @return A named list of containing (1) the sampled OD vector (named y), (2) a corresponging vector (named z)
//' the log probability density from Markov chain transition probabilities (named lq)
// [[Rcpp::export]]
List rod(IntegerVector x) {
  int S = x.size() / 2;
  IntegerVector y(((S * S) - S) / 2);
  IntegerVector z(((S * S) - S) / 2);
  NumericVector pi(S);

  int k =  0;
  int k_ijm1 = 0;
  IntegerVector kList;

  for(int j=1; j<S; ++j) {
	  // x[i] = u_i
	  // x[S+j] = v_j
    kList = IntegerVector::create();
    
    // +++ SETTING Z +++ 
	  for(int i=0; i<j; ++i) {
		  k = ij_to_id(i, j, S);
      kList.push_back(k);
		  k_ijm1= ij_to_id(i, j-1, S);

		  if(i==j-1) z(k) = x(i);
		  else z(k) = z(k_ijm1) - y(k_ijm1);
      }

    // +++ SETTING Y +++
    y[kList] = ztoy(z[kList], x(S+j));
  }

  // Markov Chain Transition Probabilities
  NumericVector w;
  NumericVector v;
  NumericVector pi_1;
  NumericVector pi_2;
  w = load(x);
  v = x[seq(S+1, x.size()-1)];
  pi_1 = -1.0 * lchoose(w[seq(0, w.size()-2)], v[seq(1,v.size()-1)]);
  pi_2 = lchoose(as<NumericVector>(z), as<NumericVector>(y));

  double lq = sum(pi_1) + sum(pi_2);

  return List::create(Named("y") = y,
		  Named("z") = z,
      Named("lq") = lq);
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
//' 
// [[Rcpp::export]]
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

      double denom = accu(exp_vals);
      for (int k = 0; k < len; ++k) {
        Lambda(start + k, col) = exp_vals(k) / denom;
      }
    }

    start += len;
  }

  return Lambda;
}

// [[Rcpp::export]]
arma::mat log_likelihood(NumericVector y, NumericVector x, NumericMatrix phi, NumericMatrix psi, double rho) {

  // perhaps best is to have this function
  // evaluate the likelihood of a single y^n
  // instead of all trips together -> needed for metropolis-hastings
  //
  // then have another function which computes the overall likelihood
  // as the product of all p(y^n | lambda^n)

  int S = x.size() / 2;
  // alternative arma::uvec (for unsigned int vectors)
  arma::vec ay = as<arma::vec>(y);
  arma::vec ax = as<arma::vec>(x);

  arma::mat aphi = as<arma::mat>(phi);
  arma::mat apsi = as<arma::mat>(psi);

  // columns -> trips
  // rows -> stops (S-2)
  arma::mat G = aphi * apsi.t();
  arma::mat lambda(G.n_rows, G.n_cols);
  
  // set lambda_S-1,S = 1
  for(arma::uword n = 0; n < G.n_cols; ++n) {
    for(int i = 0; i < S; ++i) {
      // lambda(i) = softmax(as<NumericVector>(G.row(i) * rho));
      arma::umat idx;
      idx = sub2ind(size(lambda), i_to_id(i, n, S));
      lambda.elem(idx) = softmax(G.col(n) * rho);
    }
  }

  // Eq (6.11)
  // for each bus trip
  // double log_sum_u = sum(lfactorial(x[seq(0, S-1)]));
  // double log_sum_mult = sum(log(pow(lambda, y)) - lfactorial(y));

  // arma::vec likelihood(lambda.n_cols);
  // for(arma::uword n = 0; n < lambda.n_cols; ++n) {
  //   log_sum_u = sum(lfactorial(x[seq(0, S-1)]));
  //   log_sum_mult = sum(log(pow(lambda[,n], y)) - lfactorial(y));
  // }


  // return log_sum_u + log_sum_mult;
  return lambda;
}
