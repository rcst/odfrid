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
  double gamma = randu<double>();
  double log_c = log_likelihood() + log(gamma);
  double theta = randu(distr_param(0.0, 2.0 * datum::pi));
  double theta_min = theta - (2.0 * datum::pi);
  double theta_max = theta;
  bool fQuit = false;

  // Rcpp::Rcout << "In-Psi L(): " << log_likelihood() << std::endl;
  for(uword d = 0; d < D; ++d) {
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
  }

  // update global Psi
  psi = apsi;
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
  uword D = phi.n_cols;
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
  vec nu;

  // helper - because 
  uvec d_vec;

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
    nu = mvnrnd(zeros(L), eye(size(L, L)));

    gamma = randu<double>();
    log_c = log_likelihood() + log(gamma);
    theta = randu(distr_param(0.0, 2.0 * datum::pi));
    theta_min = theta - (2.0 * datum::pi);
    theta_max = theta;

    // block indexes
    i_block = K(i, span(i+1, S-1));
    vec new_phi_d;

    // Rcpp::Rcout << ">>> Phi = " << std::endl << phi << std::endl; 

    for(uword d = 0; d < D; ++d) {
      d_vec = {d};
      trials = 0;
      fQuit = false;
      // Rcpp::Rcout << "Sampling Column " << d << std::endl;
      while (!fQuit & (trials < max_trials)) {
        new_phi_d = phi(i_block, d_vec) * std::cos(theta) + nu * std::sin(theta);
        // Rcpp::Rcout << ">>> Trial " << trials << std::endl << " Phi* (block) = " << std::endl << new_phi_d << " Phi (block) = " << std::endl << phi(i_block, d_vec) << std::endl << "<<<" << std::endl;
        aphi(i_block, d_vec) = new_phi_d;
        if(log_likelihood(d, aphi) > log_c) {
          fQuit = true;
        } else {
          // Shrink the sampling range and try a new point
          if(theta <= 0) theta_min = theta;
          else theta_max = theta;
          theta = randu(distr_param(theta_min, theta_max));
          ++trials;
        }
      }
    }

  }
  phi = aphi;
}

// // [[Rcpp::export]]
// void test(arma::uword S, arma::uword i) {
//   uword D = S * (S-1) / 2;
//   umat K(S, S);
//   uvec eids = trimatl_ind(size(K), -1);
//   K.elem(eids) = linspace<uvec>(1, D, D);
//   K = K.st();
// 
//   Rcpp::Rcout << K.rows(i, i).cols(i, i) << endl;
// 
//   Rcpp::Rcout << K << endl;
//   Rcpp::Rcout << K(i, span(i+1, S-1)) << std::endl;
// }

// // [[Rcpp::export]]
// arma::cube test() {
//   cube x(10, 10, 10, fill::randu);
//   return x;
// }
//

// [[Rcpp::export]]
arma::mat test() {
  vec M(5, fill::randu);

  mat B(5, 5, fill::randu);
  mat C = B.t() * B;

  mat X = mvnrnd(M, C, 100);
  return X;
}
