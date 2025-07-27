// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "odfrid/math.h"
#include "odfrid/likelihood.h"
#include "odfrid/od_sampler.h"
#include "odfrid/samplers.h"

#include<chrono>

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
  arma::vec L(sample);
  
  // timing
  arma::vec od_t(100);
  arma::vec phi_t(100);
  arma::vec psi_t(100);
  arma::vec rho_t(100);


  sample_od(true);

  auto start_t = std::chrono::high_resolution_clock::now();
  auto end_t = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end_t - start_t;

  arma::uword j = 0;
  for(arma::uword i = 0; i < (warmup + sample); ++i) {

    if(i >= warmup) {
      if(i % print_n == 0)
        Rcpp::Rcout << "sampling: " <<  
          i << " [" << warmup + sample << "]" << std::endl;
    } else {
      if(i % print_n == 0)
        Rcpp::Rcout << "warm-up: " <<  
          i << " [" << warmup + sample << "]" << std::endl;
    }

    // update alighting probabilities
    G = phi * psi.t();
    lbd = matrix_softmax(G, rho, S);

    // timiing estimates
    if(i < od_t.n_elem) {
      start_t = std::chrono::high_resolution_clock::now();
      sample_od_parallel();
      end_t = std::chrono::high_resolution_clock::now();
      duration = end_t - start_t;
      od_t(i) = duration.count();
    } else
      sample_od_parallel();

    if(i < phi_t.n_elem) {
      start_t = std::chrono::high_resolution_clock::now();
      ess_phi();
      end_t = std::chrono::high_resolution_clock::now();
      duration = end_t - start_t;
      phi_t(i) = duration.count();
    } else
      ess_phi();

    if(i < psi_t.n_elem) {
      start_t = std::chrono::high_resolution_clock::now();
      ess_psi(K);
      end_t = std::chrono::high_resolution_clock::now();
      duration = end_t - start_t;
      psi_t(i) = duration.count();
    } else
      ess_psi(K);

    if(i < rho_t.n_elem) {
      start_t = std::chrono::high_resolution_clock::now();
      ss_rho(0.001);
      end_t = std::chrono::high_resolution_clock::now();
      duration = end_t - start_t;
      rho_t(i) = duration.count();
    } else
      ss_rho(0.001);

    if(i == od_t.n_elem) {
      Rcpp::Rcout << std::endl << 
        "OD-sampling steps took " << arma::mean(od_t) << 
        "s on average" << std::endl;
      Rcpp::Rcout << 
       "Phi sampling steps took " << arma::mean(phi_t) << 
       "s on average" << std::endl;
      Rcpp::Rcout << 
        "Psi sampling steps took " << arma::mean(psi_t) << 
        "s on average" << std::endl;
      Rcpp::Rcout << 
        "Rho sampling steps took " << arma::mean(rho_t) << 
        "s on average" << std::endl;
      Rcpp::Rcout << 
        "Total Expected Sampling Duration is " << 
        (arma::mean(od_t) + 
         arma::mean(phi_t) + 
         arma::mean(psi_t) + 
         arma::mean(rho_t)) * (sample + warmup) / 60.0 << 
        " min." << std::endl << std::endl;
    }
    
    // Collect sampling data
    if(i >= warmup) {
      Y.slice(j) = y;
      Lambda.slice(j) = lbd;
      Phi.slice(j) = phi;
      Psi.slice(j) = psi;
      Rho(j) = rho;
      L(j) = log_likelihood();

      ++j;
    } 
  }

  return Rcpp::List::create(
      Rcpp::Named("x") = x,
      Rcpp::Named("t") = t,
      Rcpp::Named("y") = Y,
      Rcpp::Named("lambda") = Lambda,
      Rcpp::Named("phi") = Phi,
      Rcpp::Named("psi") = Psi,
      Rcpp::Named("rho") = Rho,
      Rcpp::Named("lq") = L);
}
