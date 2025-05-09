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

//' OD-vector(s) index mapping function
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

//' No. passengers on the bus immediatly after each stops
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

//' Sample OD vector y from vector of load vector z
//'
//' @param z numeric vector of no. passengers approaching a single specific stop
//' @param v double no. alighters at that specific stop
// [[Rcpp::export]]
IntegerVector ztoy(IntegerVector z, double v) {
  std::vector<double> y = uniformSimplexSample(z.length(), v);
  std::vector<int> iy = roundWithPreservedSum(y);
  iy = adjustWithCaps(iy, as<std::vector<int>>(z));
  return wrap(iy);
}

//' Conditional Sampling of OD vectors
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
  double lq = 0.0;
  w = load(x);
  v = x[seq(S+1, x.size()-1)];
  pi_1 = -1.0 * lchoose(w[seq(0, w.size()-2)], v[seq(1,v.size()-1)]);
  pi_2 = lchoose(as<NumericVector>(z), as<NumericVector>(y));

  for(const auto &x : pi_1) lq += x;
  for(const auto &x : pi_2) lq += x;

  return List::create(Named("y") = y,
		  Named("z") = z,
      Named("lq") = lq);
}

// [[Rcpp::export]]
NumericVector softmax(const NumericVector &x) {
  NumericVector ex;
  ex = exp(x);
  return ex / (1.0 * sum(ex));
}

