#include "RcppArmadillo.h"
using namespace arma;

vec log_choose_vec(const uvec& n, const uvec& k) {
  if (n.n_elem != k.n_elem) {
    Rcpp::stop("Vectors 'n' and 'k' must have the same length.");
  }

  vec n_d = conv_to<vec>::from(n);
  vec k_d = conv_to<vec>::from(k);
  vec nk_d = n_d - k_d;

  vec result(n_d.n_elem, fill::none);

  // Logical masks
  uvec is_valid = (k >= 0) && (k <= n);
  uvec is_equal = (n == k);
  uvec is_zero  = (k == 0);
  uvec compute  = is_valid && is_equal == 0 && is_zero == 0;

  // Set results
  result.fill(-datum::inf);  // default for invalid cases

  // n == k → log(choose(n, n)) = 0
  result.elem(find(is_equal)).zeros();

  // k == 0 → log(choose(n, 0)) = 0
  result.elem(find(is_zero)).zeros();

  // General case
  uvec idx = find(compute);
  if (!idx.is_empty()) {
    vec n_sub  = n_d.elem(idx);
    vec k_sub  = k_d.elem(idx);
    vec nk_sub = nk_d.elem(idx);

    vec logc = lgamma(n_sub + 1.0) - lgamma(k_sub + 1.0) - lgamma(nk_sub + 1.0);
    result.elem(idx) = logc;
  }

  return result;
}

mat log_choose_mat(const umat& N, const umat& K) {
  if (size(N) != size(K)) {
    Rcpp::stop("N and K must be the same size.");
  }

  int n_rows = N.n_rows;
  int n_cols = N.n_cols;
  mat result(n_rows, n_cols, fill::none);

  for (int col = 0; col < n_cols; ++col) {
    uvec n_col = N.col(col);
    uvec k_col = K.col(col);

    // Use the optimized log_choose_vector function
    vec logc = log_choose_vec(n_col, k_col);

    result.col(col) = logc;
  }

  return result;
}

//' Compute softmax column-wise with fixed numerator = 1 for last entry of each segment
//' 
//' NOTE: Lambda should be returned so that the entry from second to last to last is always 1
//' 
//' @param G: M x N matrix (G = Phi %*% t(Psi))
//' @param rho: temperature parameter
//' @return: M x N matrix Lambda of softmax values
mat matrix_softmax(const mat& G, double rho, uword S) {
  // G has one row less than lambda - ie, G_S-1
  uword M = G.n_rows;
  uword N = G.n_cols;

  mat Lambda(M+1, N, fill::ones);
  uword start = 0;
  uword end = 0;
  uword len = 0;

  for (uword i = 0; i < S-2; ++i) {
    // int len = segment_lengths(i);
    len = S - i - 1;
    end = start + len;

    for (uword col = 0; col < N; ++col) {
      vec g_block = G.submat(start, col, end - 1, col);
      vec scaled = rho * g_block;

      vec exp_vals(len);
      exp_vals(span(0, len - 2)) = 
        exp(scaled(span(0, len - 2)) - max(scaled));  // stabilization

      double denom = 1.0 + accu(exp_vals);
      exp_vals(len - 1) = 1.0;

      for (uword k = 0; k < len; ++k) {
        Lambda(start + k, col) = exp_vals(k) / denom;
      }
    }

    start += len;
  }

  return Lambda;
}
