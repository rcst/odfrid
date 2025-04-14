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
  // NOTE: static vars do NOT get reinitialized 
  // between subsequent function calls from R
  // ie., need main function that initializes

  // Need to use pointers to carry over the current set of parameters
  // idea: wrap everything into a struct and pass it around via a pointer
  // use a "main" function to create this struct (on the heap?)
  NumericVector y(((S * S) - S) / 2);
  NumericVector y_curr(((S * S) - S) / 2);
  NumericVector z(((S * S) - S) / 2);
  NumericVector pi(S);

  // diagnostic tracking
  IntegerVector vi(((S * S) - S) / 2);
  IntegerVector vj(((S * S) - S) / 2);
  IntegerVector vk(((S * S) - S) / 2);
  NumericVector bi(((S * S) - S) / 2);
  NumericVector aj(((S * S) - S) / 2);


  // initialize pi already here with first factor
  // pi(0) = 0 since w_0 = 0

  int k =  0;
  int yv_check = 0;
  int ytmp = 0;
  for(int j=1; j<S; ++j) {
    // do {
      yv_check = 0;
      for(int i=0; i<j; ++i) {
        k = ij_to_id(i, j, S);
	vi(k) = i;
	vj(k) = j;
	vk(k) = k;
	bi(k) = x(i);
	aj(k) = x(S + j);
        // z_j,j+1 = u_j
        // that's always the last iteration
        // needs to checking
        if(j==i+1) z(k) = x(i); 
        else z(k) = z(k-1) - y_curr(k-1);

        // alight all at the last stop
        // NOTE: Maybe faster to sample 
        // from remainder of v_j and avoid 
        // constraint check?

        // this could also be avoided
        // and done in bulk after the current
        // loop
        // if(j == S-1) y_curr(k) = z(k);
        // else y_curr(k) = (int)runif(1, 0, z(k)+1)[0];

        // yv_check += y_curr(k);
        // check for early exit
        // still bad, because it loses all
        // previous iterations
        // better to loop over random draw:
        if(j == S-1) 
          y_curr(k) = z(k);
        else if(z(k) == 0)
		y(k) = 0;
	else {
          do {
            ytmp = (int)runif(1, 0, z(k)+1)[0];
          } while ((yv_check + ytmp) > x(S+j));
          y_curr(k) = ytmp;
        }
        yv_check += ytmp;

        // if(yv_check > x(S+j)) {
        //   break;
        // }
	// Rcout << "i: " << i << " j: " << j << " k: " << k << " z: " << z(k) << " y: " << y_curr(k) << std::endl;
      }
    // } while(yv_check != x(S+j));
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
  pi2 = choose(z, y_curr);
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

  // accept or not?
  y = y_curr;

  // diagnostic data.frame
  DataFrame dg = DataFrame::create(
		  Named("i") = vi, 
		  Named("j") = vj, 
		  Named("k") = vk, 
		  Named("z") = z, 
		  Named("y") = y, 
		  Named("ui") = bi, 
		  Named("vj") = aj);

  return List::create(Named("y") = y,
                      Named("z") = z,
                      Named("pi") = pi,
                      Named("q") = q,
		      Named("dg") = dg);
}

// [[Rcpp::export]]
NumericVector softmax(const NumericVector &x) {
  NumericVector ex;
  ex = exp(x);
  return ex / (1.0 * sum(ex));
}

