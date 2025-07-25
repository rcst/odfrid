#include "RcppArmadillo.h"
#include "odfrid/globals.h"
#include "odfrid/math.h"
#include "odfrid/likelihood.h"

using namespace arma;

mat departure_times_covariance_matrix(const vec& t, double sigma, double l) {
  mat T1 = repmat(t, 1, t.n_elem);
  mat T2 = repmat(t.t(), t.n_elem, 1);
  mat K = std::pow(sigma, 2) * exp(square(T1 - T2) / (-2.0 * std::pow(l, 2)));      

  return K;
}

//' ess_psi - Elliptical Slice Sampling of Psi
//'
//' One call to this function wll update the global matrix Psi based on the 
//' current (global) set of parameters.
//'
//' @param k covariance matrix for column psi_d
void ess_psi(mat& K) {
  // uword N = psi.n_cols;
  // uword D = psi.n_rows;
  uword N = psi.n_rows;
  uword D = psi.n_cols;
  vec nu = mvnrnd(zeros(N), K);
  // copy
  mat apsi = psi;
  double gamma;
  double log_c;
  double theta;
  double theta_min;
  double theta_max;
  bool fQuit = false;

  // This could be accelerated by boardcasting i.e. .each_col()
  // Qestion to Xiaxio: Is the likelihood evaluated based on old columns
  // ... or already updated previous columns?
  for(uword d = 0; d < D; ++d) {
    gamma = randu<double>();
    log_c = log_likelihood() + log(gamma);
    theta = randu(distr_param(0.0, 2.0 * datum::pi));
    theta_min = theta - (2.0 * datum::pi);
    theta_max = theta;
    fQuit = false;
    while (!fQuit) {
      vec new_psi_d = psi.col(d) * std::cos(theta) + nu * std::sin(theta);
      apsi.col(d) = new_psi_d;
      if(log_likelihood(apsi, d) <= log_c) {
        // Shrink the sampling range and try a new point
        if(theta <= 0) theta_min = theta;
        else theta_max = theta;
        theta = randu(distr_param(theta_min, theta_max));
      } else {
        fQuit = true;
      }
    }
  // update global Psi
  psi = apsi;
  }
}

void ss_rho(double eps) {
  double gamma = randu<double>();
  double log_p_rho = log_normpdf(rho, log(0.1), 1.0);
  double log_c = log_likelihood() + log_p_rho + log(gamma);

  double kappa = randu<double>(distr_param(0.0, eps));
  double rho_min = rho - kappa;
  double rho_max = rho_min + eps;

  bool fQuit = false;
  double new_rho  = 0.0;

  while(!fQuit) {
    new_rho = randu<double>(distr_param(rho_min, rho_max));

    if(log_likelihood(new_rho) + log_normpdf(new_rho, log(0.1), 1.0) > log_c) {
      fQuit = true; 
    } else {
      if(new_rho < rho) rho_min = new_rho;
      else rho_max = new_rho;
    }
  }

  rho = new_rho;
}

void ess_phi() {
  // one sampling steps updates whole matrix
  // updates in a block-wise
  // sample each column d for an i
  // block is phi_i for all d
  // instead of column-wise
  // uword D = phi.n_cols;
  uword N = phi.n_rows; 
  uword S = od_size_to_nstops(N+1);
  uword M = S * (S - 1) / 2;
  uword L = 0;

  // copy
  mat aphi = phi;

  double gamma = 0;
  double log_c = 0;
  double theta = 0;
  double theta_min = 0;
  double theta_max = 0;
  bool fQuit = false;
  mat nu;

  // helper - because 
  //uvec d_vec;

  // tracking indices of vectorized matrices 
  umat K(S, S);
  Row<uword> i_block;
  uvec eids = trimatl_ind(size(K), -1);
  K.elem(eids) = linspace<uvec>(0, M-1, M);
  K = K.st();

  uword max_trials = 100;
  uword trials = 0;

  for(uword i = 0; i < (S-2); ++i) {
    L = S-i-1;
    nu = mvnrnd(zeros(L), eye(size(L, L)), phi.n_cols);

    gamma = randu<double>();
    log_c = log_likelihood() + log(gamma);
    theta = randu(distr_param(0.0, 2.0 * datum::pi));
    theta_min = theta - (2.0 * datum::pi);
    theta_max = theta;

    // block indexes
    i_block = K(i, span(i+1, S-1));
    mat new_phi_d;

    trials = 0;
    fQuit = false;
    while (!fQuit & (trials < max_trials)) {
      new_phi_d = phi.rows(i_block) * std::cos(theta) + nu * std::sin(theta);
      
      // replace several blocks
      aphi.rows(i_block) = new_phi_d;
      if(log_likelihood(1, aphi) > log_c) {
        fQuit = true;
      } else {
        // Shrink the sampling range and try a new point
        if(theta <= 0) theta_min = theta;
        else theta_max = theta;
        theta = randu(distr_param(theta_min, theta_max));
        ++trials;
        if(trials == max_trials)
          Rcpp::Rcout << "<<< Max. no. trials reached! >>>" << std::endl;
      }
    }
    phi.rows(i_block) = aphi.rows(i_block);
    }
  }

// [[Rcpp::export]]
void test() {
  umat y(300, 2, fill::value(3));
  mat lbd(300, 100, fill::randu);

  vec lq = sum((y % log(lbd)) - lgamma(y));

  Rcpp::Rcout << lq << std::endl;
}
