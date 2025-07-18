#ifndef ODMATH_H 
#define ODMATH_H

arma::vec log_choose_vec(const arma::uvec& n, const arma::uvec& k);
arma::mat log_choose_mat(const arma::umat& N, const arma::umat& K);
arma::vec softmax(const arma::vec &x);
arma::mat matrix_softmax(const arma::mat& G, double rho);
arma::uword od_size_to_nstops(arma::uword N);

#endif

