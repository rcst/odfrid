#ifndef SAMPLER_H
#define SAMPLER_H

arma::imat load(arma::imat& x);
Rcpp::List rod(arma::imat& x);
arma::umat odform2sub(arma::umat& K, int j, int c);
arma::ivec ztoy(arma::ivec z, double v);
arma::mat departure_times_covariance_matrix(const arma::vec& t, double sigma, double l);
arma::mat ess_psi(arma::mat& K, arma::uword& d);
double ss_rho(double eps);
arma::mat ess_phi(arma::uword& d);

#endif

