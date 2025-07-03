#ifndef ODMATH_H 
#define ODMATH_H

arma::vec log_choose_vec(const arma::ivec& n, const arma::ivec& k);
arma::mat log_choose_mat(const arma::imat& N, const arma::imat& K);
arma::vec softmax(const arma::vec &x);
arma::mat matrix_softmax(const arma::mat& G, double rho);

#endif

