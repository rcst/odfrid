#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

arma::vec log_likelihood(const arma::imat& y, const arma::imat& x, const arma::mat& phi, const arma::mat& psi);
arma::vec log_likelihood(const arma::imat& y, const arma::imat& x);
arma::vec log_likelihood(const arma::mat& psi);

#endif

