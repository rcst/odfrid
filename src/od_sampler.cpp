#include "RcppArmadillo.h"
#include "odfrid/globals.h"
#include "odfrid/math.h"
#include "odfrid/uniform_simplex.h"
#include "odfrid/likelihood.h"

using namespace arma;

//' load - No. passengers on the bus immediatly after each stops
//'
//' @param x column vector of boardings and alightings
//' @return vector of passengers loadings immediatly after each stops
// [[Rcpp::export]]
arma::umat load(arma::umat& x) {
  uword S = x.n_rows / 2.0;
  umat w(S, x.n_cols);

  for(uword c = 0; c < x.n_cols; ++c) {
    for(uword i = 0; i<S; ++i) {
      if(i==0) w(i, c) = x(i, c) - x(S+i, c);
      else w(i, c) = w(i-1, c) + x(i, c) - x(S+i, c);
    }
  }

  return w;
}

//' ztoy - capped uniform simplex sampling 
//'
//' Helper function that generates an integer vector with constant total sum and varying caps
//'
//' @param z numeric vector of no. passengers approaching a single specific stop
//' @param v double no. alighters at that specific stop
//' @return An integer vector of same length as z whos values are all smaller-or-equal to z and it's sum is equal to v
uvec ztoy(uvec z, double v) {
  std::vector<double> y = uniformSimplexSample(z.n_elem, v);
  std::vector<unsigned int> uy = roundWithPreservedSum(y);
  // iy = adjustWithCaps(iy, as<std::vector<int>>(z));
  //Rcpp::Rcout << "uy = " << std::endl << uvec(uy) << std::endl << "z = " << std::endl << z << "v = " << v << std::endl;
  uy = adjustWithCaps(uy, conv_to<std::vector<unsigned int>>::from(z));
  uvec r(uy.data(), uy.size(), true);
  return r;
}

umat odform2sub(umat& K, int j, int c) {
  uvec C(j);
  C.fill(c);

  return join_rows(K(span(0, j-1), j), C).st();
}

uvec odform2sub(umat& K, int j) {
  return K(span(0, j-1), j);
}

//' rod - Conditional Sampling of OD vectors
void sample_od(bool force_accept = false) {
  uword N = x.n_cols;
  uword S = x.n_rows / 2;
  uword M = S * (S-1) / 2;

  umat u = x.rows(0, S-1);
  umat v = x.rows(S, x.n_rows - 1);
  // y, z and q candidate
  uvec yc(M);
  uvec zc(M);
  double qc;
  double p_acc = 0.0;

  // Markov Chain Transition Probabilities
  umat w = load(x);
  // "constant" factor from data
  // no need to re-compute everytime
  const mat pi_1 = -1.0 * log_choose_mat(w.rows(0, w.n_rows - 2), v.rows(1, v.n_rows - 1));

  // tracking indices of vectorized matrices 
  umat K(S, S);
  uvec eids = trimatl_ind(size(K), -1);
  K.elem(eids) = linspace<uvec>(0, M-1, M);
  K = K.st();

  uvec klist;
  uword k =  0;
  uword k_ijm1 = 0;

  for(uword c = 0; c < N; ++c) {
    for(uword j = 1; j < S; ++j) {
      // +++ SETTING Z +++ 
      for(uword i = 0; i < j; ++i) {
        k = K(i, j);
        k_ijm1 = K(i, j-1);
        if(i==j-1) zc(k) = u(i, c);
        else zc(k) = zc(k_ijm1) - yc(k_ijm1);
      }

      // +++ SETTING Y +++
      klist = odform2sub(K, j);
      //Rcpp::Rcout << ">>> j = " << j << " c = " << c << std::endl;
      yc.elem(klist) = ztoy(zc.elem(klist), v(j, c));
    }
    // decide whether to accept or dismiss candidate y, z and q
    qc = sum(pi_1.col(c)) + sum(log_choose_vec(zc, yc));
    // p_acc = (accu(log_likelihood(yc)) * q(c)) / (accu(log_likelihood()) * qc);

    if(force_accept) {
      y.col(c) = yc;
      z.col(c) = zc;
      q(c) = qc;
    } else {
      p_acc = (log_likelihood(yc, x(span::all, c)) + q(c)) - 
        (log_likelihood(y(span::all, c), x(span::all, c)) + qc);
      if(p_acc >= 0.0) {
        // accept right away
        y.col(c) = yc;
        z.col(c) = zc;
        q(c) = qc;
      } else {
        // accept with random probability
        // comparing on log-scale
        if(log(randu()) < p_acc) {
          y.col(c) = yc;
          z.col(c) = zc;
          q(c) = qc;
        }
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::List rod(const arma::umat& x) {
  uword N = x.n_cols;
  uword S = x.n_rows / 2;
  uword M = S * (S-1) / 2;

  umat u = x.rows(0, S-1);
  umat v = x.rows(S, x.n_rows - 1);
  // y, z and q candidate
  uvec yc(M);
  uvec zc(M);
  // double qc;
  // double p_acc = 0.0;

  //result
  umat y(M, N);
  umat z(M, N);

  // Markov Chain Transition Probabilities
  // umat w = load(x);
  // "constant" factor from data
  // no need to re-compute everytime
  // const mat pi_1 = -1.0 * log_choose_mat(w.rows(0, w.n_rows - 2), v.rows(1, v.n_rows - 1));

  // tracking indices of vectorized matrices 
  umat K(S, S);
  uvec eids = trimatl_ind(size(K), -1);
  K.elem(eids) = linspace<uvec>(0, M-1, M);
  K = K.st();

  uvec klist;
  uword k =  0;
  uword k_ijm1 = 0;

  for(uword c = 0; c < N; ++c) {
    for(uword j = 1; j < S; ++j) {
      // +++ SETTING Z +++ 
      for(uword i = 0; i < j; ++i) {
        k = K(i, j);
        k_ijm1 = K(i, j-1);
        if(i==j-1) zc(k) = u(i, c);
        else zc(k) = zc(k_ijm1) - yc(k_ijm1);
      }

      // +++ SETTING Y +++
      klist = odform2sub(K, j);
      yc.elem(klist) = ztoy(zc.elem(klist), v(j, c));
    }
      y.col(c) = yc;
      z.col(c) = zc;
  }

  return(Rcpp::List::create(
        Rcpp::Named("y") = y,
        Rcpp::Named("z") = z));
}
