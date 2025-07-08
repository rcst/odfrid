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

//' @param k covariance matrix for column psi_d
//' @param psi mapping-factor matrix
//' @param d the index of column of matrix psi to be sampled
mat ess_psi(mat& K, uword& d) {
  uword N = psi.n_cols;
  vec nu = mvnrnd(zeros(N), K);

  // copy
  mat apsi = psi;

  double gamma = randu<double>();

  double log_c = accu(log_likelihood()) + log(gamma);
  double theta = randu(distr_param(0.0, 2.0 * datum::pi));
  double theta_min = theta - (2.0 * datum::pi);
  double theta_max = theta;


  bool fQuit = false;
  while (!fQuit) {
    vec new_psi_d = psi.col(d) * std::cos(theta) + nu * std::sin(theta);
    apsi.col(d) = new_psi_d;
    if(accu(log_likelihood(apsi, d)) <= log_c) {
      // Shrink the sampling range and try a new point
      if(theta <= 0) theta_min = theta;
      else theta_max = theta;

      theta = randu(distr_param(theta_min, theta_max));
    } else fQuit = true;
  }

  return apsi;
}

double ss_rho(double eps) {
  double gamma = randu<double>();
  double log_p_rho = log_normpdf(rho, log(0.1), 1.0);
  double log_c = accu(log_likelihood()) + log_p_rho + log(gamma);

  double kappa = randu<double>(distr_param(0.0, eps));
  double rho_min = rho - kappa;
  double rho_max = rho_min + eps;

  bool fQuit = false;
  double new_rho  = 0.0;

  while(!fQuit) {
    new_rho = randu<double>(distr_param(rho_min, rho_max));

    if(accu(log_likelihood(new_rho)) + log_normpdf(new_rho, log(0.1), 1.0) > log_c) {
      fQuit = true; 
    } else {
      if(new_rho < rho) rho_min = new_rho;
      else rho_max = new_rho;
    }
  }

  return new_rho;
}

mat ess_phi(uword& d) {
  // update block-wise
  // sample each column d for a i
  // block is phi_i for all d
  // instead of column-wise
  uword N = phi.n_cols;
  vec nu = mvnrnd(zeros(N), eye(size(N, N)));

  // copy
  mat aphi = phi;

  double gamma = randu<double>();

  double log_c = accu(log_likelihood()) + log(gamma);
  double theta = randu(distr_param(0.0, 2.0 * datum::pi));
  double theta_min = theta - (2.0 * datum::pi);
  double theta_max = theta;


  bool fQuit = false;
  while (!fQuit) {
    vec new_phi_d = phi.col(d) * std::cos(theta) + nu * std::sin(theta);
    aphi.col(d) = new_phi_d;
    if(accu(log_likelihood(aphi, d)) <= log_c) {
      // Shrink the sampling range and try a new point
      if(theta <= 0) theta_min = theta;
      else theta_max = theta;

      theta = randu(distr_param(theta_min, theta_max));
    } else fQuit = true;
  }

  return aphi;
}
