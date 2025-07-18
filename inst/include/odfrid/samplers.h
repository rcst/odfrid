#ifndef SAMPLER_H
#define SAMPLER_H

arma::mat departure_times_covariance_matrix(const arma::vec& t, double sigma, double l);
void ess_psi(arma::mat& K);
double ss_rho(double eps);
void ess_phi();
void test(arma::uword S, arma::uword i);

#endif

