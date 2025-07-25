#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

double log_likelihood(const arma::umat& y, const arma::umat& x, const arma::mat& phi, const arma::mat& psi, const double& rho, const arma::umat& A);

double log_likelihood(const arma::uvec& y, const arma::uvec& x, const arma::uword n);
double log_likelihood(const arma::umat& y, const arma::umat& x);
double log_likelihood(const arma::mat& psi, const arma::uword d);
double log_likelihood(const arma::uword d, const arma::mat& psi);
double log_likelihood(const double& rho);
double log_likelihood();

#endif

