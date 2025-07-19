// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "odfrid/math.h"
#include "odfrid/likelihood.h"
#include "odfrid/od_sampler.h"
#include "odfrid/samplers.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]

// define global model parameters
arma::umat x;
arma::vec t;
arma::umat A;
arma::umat y;
arma::umat z;
arma::vec q;
arma::mat psi;
arma::mat phi;
double rho;


//' routing_matrix - construct routing matrix from no. stops (S)
//'
//' @param S length-one integer denoting the number of stops
//' @return An (armadillo) integer matrix
// [[Rcpp::export]]
arma::umat routing_matrix(arma::uword s) {

  arma::uword m = s * (s - 1) / 2;
  arma::umat A(2*s, m);

  // column index variable of A
  arma::uword j = 0;

  for(arma::uword i = 0; i < s; ++i) {
    for(arma::uword k = i + 1; k < s; ++k) {
      j = s * i - (i + 1) * (i + 2) / 2 + k;
      A(i, j) = 1.0;
    }
  }

  for(arma::uword i = s; i < 2*s; ++i) {
    for(arma::uword k = 0; k < (i + 1) - s - 1; ++k) {
      j = s * (k - 1) - (k + 1) * (k + 2) / 2 + i;
      A(i, j) = 1.0;
    }
  }

  return A;
}

// [[Rcpp::export]]
Rcpp::List model_sample(
    const arma::umat& ax, 
    const arma::vec& dep_time,
    const arma::uword sample, 
    const arma::uword warmup, 
    const arma::uword D,
    const arma::uword print_n = 100) {

  const arma::uword S = ax.n_rows / 2;
  const arma::umat& u = ax(arma::span(0, S-1), arma::span::all); 
  const arma::umat& v = ax(arma::span(S, (2*S)-1), arma::span::all);

  const arma::uword N = u.n_cols;
  const arma::uword M = S * (S - 1) / 2;

  // allocate global variables
  A = routing_matrix(S);
  y = arma::umat(M, N);
  z = arma::umat(M, N);
  q = arma::vec(N);

  t = dep_time;
  x = join_cols(u, v);
  arma::mat K = departure_times_covariance_matrix(t, 1, 3600.0);

  rho = 1.0;
  phi = arma::mvnrnd(arma::zeros(M-1), arma::eye(arma::size(M-1, M-1)), D);
  psi = arma::mvnrnd(arma::zeros(N), arma::eye(arma::size(N, N)), D);


  arma::mat G;
  arma::mat lbd;

  // output
  arma::ucube Y(M, N, sample);
  arma::cube Lambda(M, N, sample);
  arma::cube Phi(phi.n_rows, phi.n_cols, sample);
  arma::cube Psi(psi.n_rows, psi.n_cols, sample);
  arma::vec Rho(sample);
  


  sample_od(true);

  arma::uword j = 0;
  for(arma::uword i = 0; i < (warmup + sample); ++i) {
    // update alighting probabilities
    G = phi * psi.t();
    lbd = matrix_softmax(G, rho);

    sample_od();

    if(i >= warmup) {
      if(i % print_n == 0)
        Rcpp::Rcout << "sampling: " <<  i << " [" << warmup + sample << "]" << std::endl;
      
      // collect OD vector and alighting probabilities
      Y.slice(j) = y;
      Lambda.slice(j) = lbd;
      Phi.slice(j) = phi;
      Psi.slice(j) = psi;
      Rho(j) = rho;

      ++j;
    } else {
      if(i % print_n == 0)
        Rcpp::Rcout << "warm-up: " <<  i << " [" << warmup + sample << "]" << std::endl;
    }

    // Rcpp::Rcout << "Phi" << std::endl;
    ess_phi();

    // Rcpp::Rcout << "Psi" << std::endl;
    ess_psi(K);

    // Rcpp::Rcout << "Rho" << std::endl;
    ss_rho(0.001);

  }

  return Rcpp::List::create(
      Rcpp::Named("y") = Y,
      Rcpp::Named("lambda") = Lambda,
      Rcpp::Named("phi") = Phi,
      Rcpp::Named("psi") = Psi,
      Rcpp::Named("rho") = rho,
      Rcpp::Named("lq") = log_likelihood());
}
