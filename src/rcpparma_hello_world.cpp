// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

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
NumericVector load(NumericVector x) {
  int S = x.size() / 2.0;
  // std::cout << "S: " << S << std::endl;
  NumericVector w(S);

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

//' Conditional Sampling of OD vectors
//' 
//' @param x a column vector containing boarding and alighting counts of 1 bus
//' journey
//' @return A named list of containing (1) the sampled OD vector, (2) a corresponging vector (z)
//' (3) a vector of corresponding Markov chain transition probabilities
// [[Rcpp::export]]
List rod(NumericVector x) {
  // TODO when we add the acceptance step we need to add the model parameters
  // $Theta$ need to be added as a parameter
  // if running safely, replace access-operator () by [] which avoids bound check and
  // is faster
  
  int S = x.size() / 2;
  NumericVector w;
  NumericVector y(((S * S) - S) / 2);
  NumericVector z(((S * S) - S) / 2);
  NumericVector pi(S);

  int k =  0;
  int kjm1 = 0;
  int yv_check = 0;
  for(int j=1; j<S; ++j) {
      yv_check = 0;
      for(int i=0; i<j; ++i) {
        k = ij_to_id(i, j, S);
	kjm1 = ij_to_id(i, j-1, S);

	// +++ SETTING Z +++ 
	// z_j,j+1 = u_j
        if(i==j-1) z(k) = x(i);
        else z(k) = z(kjm1) - y(kjm1);

	Rcout << z(k) << "(" << k << ")" << std::endl;

	if(z(k) < 0) {
		Rcout << "neg. z!" << std::endl; 
		Rcout << "i: " << i << " j: " << j << " k: " << k << " z:" << z(k) << " zprev: " << z(kjm1) << " y: " << y(k) << " y_prev: " << y(kjm1) << " sum: " << yv_check << " ui: "<< x(j) << " vj: " << x(S+j) << std::endl;
	}

	// +++ SETTING Y +++
	// alight everyone at last stop
	if(j == S-1) y(k) = z(k);
	// if none on board, none to alight
	else if(z(k) == 0) y(k) = 0;
	// last i of given j -> alight remainder
	else if(i==j-1)	y(k) = x(S+j) - yv_check;
	else
		// sample randomly
		do y(k) = (int)runif(1, 0, z(k)+1)[0];
		while ((yv_check + y(k)) > x(S+j));


	if(y(k) < 0) {
		Rcout << "neg. y !" << std::endl; 
		Rcout << "i: " << i << " j: " << j << " k: " << k << " z:" << z(k) << " y: " << y(k) << "sum: " << yv_check << " ui: "<< x(j) << " vj: " << x(S+j) << std::endl;
	}
	yv_check += y(k);
      }
  }

  // Markov Chain Transition Probabilities
  w = load(x);
  NumericVector v;
  v = x[seq(S+1, x.size()-1)];
  pi = 1 / choose(w[seq(0, w.size()-2)], v[seq(1,v.size()-1)]);
  // using log-choose for computational efficiency
  // log_pi = -1 * lchoose(w[seq(0, w.size()-2)], v[seq(1,v.size()-1)]);
  
  // perhaps better to use lchoose
  // and product sums
  NumericVector pi2;
  pi2 = choose(z, y);
  // pi2 = lchoose(z, y);
  double pi3 = 1.0;
  // double pi3 = 0.0;
  for(int i=0; i<pi2.size(); ++i) {
    pi3 *= pi2[i];
    // pi3 += pi2[i];
  }

  double q = 1.0;
  // double q = 0.0;
  for(int i=0; i<pi.size(); ++i) {
    q *= pi[i];
    // q += pi[i];
  }

  q *= pi3; 

  return List::create(Named("y") = y,
                      Named("z") = z,
                      Named("pi") = pi,
                      Named("q") = q);
}

// [[Rcpp::export]]
NumericVector softmax(const NumericVector &x) {
  NumericVector ex;
  ex = exp(x);
  return ex / (1.0 * sum(ex));
}

