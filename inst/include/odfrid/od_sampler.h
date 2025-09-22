#ifndef OD_SAMPLER_H
#define OD_SAMPLER_H

#include "odfrid/globals.h"

arma::imat load(arma::imat& x);
arma::ivec ztoy(arma::ivec z, double v);
arma::umat odform2sub(arma::umat& K, int j, int c);
arma::uvec odform2sub(arma::umat& K, int j);
void sample_od(bool force_accept = false, arma::uword trips_begin = 0, arma::uword trips_end = y.n_cols);
void sample_od_parallel();

#endif
