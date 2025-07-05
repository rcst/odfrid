#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

arma::vec log_likelihood(const arma::imat& y, const arma::imat& x, const arma::mat& phi, const arma::mat& psi, const double& rho);
arma::vec log_likelihood(const arma::imat& y, const arma::imat& x);

// to isolate the signatures, we use dummy arguments and change the
// order of arguments
arma::vec log_likelihood(const arma::mat& psi, const arma::uword d);
arma::vec log_likelihood(const arma::uword d, const arma::mat& psi);

arma::vec log_likelihood(const double& rho);
arma::vec log_likelihood();

#endif

