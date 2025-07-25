// [[Rcpp::depends(RcppParallel)]]
#include "RcppParallel.h"
#include "RcppArmadillo.h"
#include "odfrid/globals.h"
#include "odfrid/math.h"
#include "odfrid/uniform_simplex.h"
#include "odfrid/likelihood.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

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
      p_acc = (log_likelihood(yc, x.col(c), c) + q(c)) - 
        (log_likelihood(y.col(c), x.col(c), c) + qc);
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

// Assign arma::uvec `src` to column `col_index` of RMatrix<uint32_t> `dest`
void assign_uvec_to_column(RMatrix<int>& dest, const arma::uvec& src, std::size_t col_index) {
  std::size_t n_rows = dest.nrow();

  if (src.n_elem != n_rows) {
    stop("Size mismatch: source vector and destination matrix column have different lengths.");
  }

  for (std::size_t i = 0; i < n_rows; ++i) {
    dest(i, col_index) = src[i];
  }
}

struct ODSamplingWorker: public Worker {
  // const RMatrix<int> input;
  RMatrix<int> output;

  uword N;
  uword S;
  uword M;

  umat u;
  umat v;
 
  // y, z and q candidate
  uvec yc;
  uvec zc;
  double qc;
  double p_acc;

  // masks globals
  umat z;
  uvec q;

  // Markov Chain Transition Probabilities
  umat w;
  // "constant" factor from data
  // no need to re-compute everytime
  mat pi_1;

  // tracking indices of vectorized matrices 
  umat K;

  uvec klist;

  ODSamplingWorker(IntegerMatrix& aux_y)
    : output(aux_y) {
        N = x.n_cols;
        S = x.n_rows / 2;
        // N = aux_x.ncol();
        // S = aux_x.nrow() / 2;
        M = S * (S-1) / 2;

        // zero-copy
        // x = umat(reinterpret_cast<uint32_t*>(aux_x.begin()),
        //       aux_x.nrow(),
        //       aux_x.ncol());
        // x = as<arma::umat>(aux_x);

        u = x.rows(0, S-1);
        v = x.rows(S, x.n_rows - 1);

        // y, z and q candidate
        yc = uvec(M);
        zc = uvec(M);
        qc = 0.0;
        p_acc = 0.0;

        // Markov Chain Transition Probabilities
        w = load(x);
        // "constant" factor from data
        // no need to re-compute everytime
        mat pi_1 = -1.0 * log_choose_mat(w.rows(0, w.n_rows - 2), v.rows(1, v.n_rows - 1));

        // tracking indices of vectorized matrices 
        K = umat(S, S);
        uvec eids = trimatl_ind(size(K), -1);
        K.elem(eids) = linspace<uvec>(0, M-1, M);
        K = K.st();
      }

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t c = begin; c < end-1; ++c) {
      for(uword j = 1; j < S; ++j) {
        //
        // +++ SETTING Z +++ 
        //
        for(uword i = 0; i < j; ++i) {
          if(i==j-1) zc(K(i,j)) = u(i, c);
          else zc(K(i,j)) = zc(K(i,j-1)) - yc(K(i,j-1));
        }

        //
        // +++ SETTING Y +++
        //
        klist = odform2sub(K, j);
        yc.elem(klist) = ztoy(zc.elem(klist), v(j, c));
      }

      // accept or dismiss candidate
      qc = sum(pi_1.col(c)) + sum(log_choose_vec(zc, yc));
      p_acc = (log_likelihood(yc, x.col(c)) + q(c)) - 
        (log_likelihood(y.col(c), x.col(c), c) + qc);

      if(p_acc >= 0.0) {
        // y.col(c) =  yc;
        assign_uvec_to_column(output, yc, c);
        z.col(c) = zc;
        q(c) = qc;
      } else {
        if(log(randu()) < p_acc) {
          // y.col(c) = yc;
          assign_uvec_to_column(output, yc, c);
          z.col(c) = zc;
          q(c) = qc;
        }
      }
    }
  }
};

void sample_od_parallel() {
  IntegerMatrix aux_y = wrap(y);

  ODSamplingWorker ods(aux_y);
  parallelFor(0, y.n_cols, ods);

  y = Rcpp::as<arma::umat>(aux_y);
}
