#ifndef OD_SAMPLER_H
#define OD_SAMPLER_H

arma::imat load(arma::imat& x);
arma::ivec ztoy(arma::ivec z, double v);
arma::umat odform2sub(arma::umat& K, int j, int c);
Rcpp::List rod(arma::imat& x);

#endif
